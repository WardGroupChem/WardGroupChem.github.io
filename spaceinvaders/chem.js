/* ============================================================
   chem.js — SMILES → 2D skeletal structure (SVG)

   One engine for every structure in the package, so drawing is
   consistent by construction rather than hand-tuned per molecule.
   It follows the same conventions as RDKit / ChemDraw depictions:
     • one fixed bond length for every bond in every molecule
     • 120° zigzag chains; branches at 120° (trigonal) or splayed
       evenly for quaternary centres
     • sp centres (triple bonds, cumulated C=C=C) drawn linear
     • E/Z geometry taken from / and \ in the SMILES and enforced
     • rings as regular polygons, aromatics drawn Kekulé
     • double bonds offset to the inside of a ring or the substituted
       side; centred when one end is terminal (C=O, =CH2)
     • heteroatoms labelled with their implicit H (OH, NH₂, HO–) and
       real subscripts; bonds clipped so they meet the heteroatom
       symbol itself, never the H
     • skeletal: carbons are unlabelled vertices

   Usage:  Chem.svg('C/C=C/C')              -> '<svg …>'
           Chem.pair('CC(Br)CC', 'C/C=C/C') -> [svgA, svgB] (equal heights)
   Colours use currentColor, so the structure inherits the text colour
   of whatever it sits in (works on light or dark themes).
   ============================================================ */

var Chem = (function () {
  'use strict';

  var CFG = {
    bond: 30,          // px per bond — identical everywhere
    font: 13,          // atom label size, px
    lineWidth: 1.7,
    gap: 0.17,         // spacing between lines of a multiple bond (× bond)
    trim: 0.13,        // shortening of an inner (offset) double-bond line at each end (× bond)
    margin: 5,
    fontFamily: 'Helvetica, Arial, sans-serif'
  };

  var VALENCE = { B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6], F: [1], Cl: [1], Br: [1], I: [1] };
  var DEG = Math.PI / 180;

  /* ---------------- SMILES parser ---------------- */

  function parse(smiles) {
    var atoms = [], bonds = [];
    var stack = [], prev = -1, pend = null, rings = {};
    var s = smiles, i = 0;

    function addBond(a, b, sym) {
      var arom = (sym === null || sym === ':') && atoms[a].arom && atoms[b].arom;
      bonds.push({
        a: a, b: b,
        order: sym === '=' ? 2 : sym === '#' ? 3 : 1,
        arom: arom,
        dir: sym === '/' ? 1 : sym === '\\' ? -1 : 0
      });
    }

    while (i < s.length) {
      var c = s[i];
      if (c === '(') { stack.push(prev); i++; continue; }
      if (c === ')') { prev = stack.pop(); i++; continue; }
      if ('-=#:/\\'.indexOf(c) >= 0) { pend = c; i++; continue; }
      if (c === '.') { prev = -1; i++; continue; }
      if (/[0-9%]/.test(c)) {
        var num;
        if (c === '%') { num = s.substr(i + 1, 2); i += 3; } else { num = c; i++; }
        if (rings[num]) {
          var r = rings[num];
          var sym = pend || r.bond;
          if (sym === '/' || sym === '\\') sym = null; // stereo on ring closures not supported
          addBond(r.atom, prev, sym);
          delete rings[num];
        } else {
          rings[num] = { atom: prev, bond: pend };
        }
        pend = null;
        continue;
      }

      var atom = { el: null, arom: false, bracket: false, hExp: 0, charge: 0 };
      if (c === '[') {
        var j = s.indexOf(']', i);
        var inner = s.slice(i + 1, j);
        var m = /^(\d*)([A-Z][a-z]?|[a-z][a-z]?)(@*)(H\d*)?([+-]+\d*)?/.exec(inner);
        if (!m) throw new Error('Bad bracket atom [' + inner + '] in ' + smiles);
        var sym2 = m[2];
        atom.arom = sym2[0] === sym2[0].toLowerCase();
        atom.el = sym2[0].toUpperCase() + sym2.slice(1);
        atom.bracket = true;
        atom.hExp = m[4] ? (m[4].length > 1 ? parseInt(m[4].slice(1), 10) : 1) : 0;
        if (m[5]) {
          var sign = m[5][0] === '+' ? 1 : -1;
          var digits = m[5].replace(/[+-]/g, '');
          atom.charge = sign * (digits ? parseInt(digits, 10) : m[5].length);
        }
        i = j + 1;
      } else if (s.substr(i, 2) === 'Cl' || s.substr(i, 2) === 'Br') {
        atom.el = s.substr(i, 2); i += 2;
      } else if (/[BCNOPSFI]/.test(c)) {
        atom.el = c; i++;
      } else if (/[bcnops]/.test(c)) {
        atom.el = c.toUpperCase(); atom.arom = true; i++;
      } else {
        throw new Error('Unexpected "' + c + '" in SMILES ' + smiles);
      }
      atoms.push(atom);
      var idx = atoms.length - 1;
      if (prev >= 0) addBond(prev, idx, pend);
      pend = null;
      prev = idx;
    }

    var mol = { atoms: atoms, bonds: bonds, smiles: smiles };
    buildAdj(mol);
    kekulize(mol);
    implicitH(mol);
    findRings(mol);
    return mol;
  }

  function buildAdj(mol) {
    mol.adj = mol.atoms.map(function () { return []; });
    mol.bonds.forEach(function (b, k) {
      mol.adj[b.a].push({ nb: b.b, bond: k });
      mol.adj[b.b].push({ nb: b.a, bond: k });
    });
  }

  function bondBetween(mol, a, b) {
    var list = mol.adj[a];
    for (var k = 0; k < list.length; k++) if (list[k].nb === b) return mol.bonds[list[k].bond];
    return null;
  }

  /* Aromatic lowercase SMILES → alternating single/double (Kekulé) by
     backtracking perfect matching over the aromatic bonds. */
  function kekulize(mol) {
    var atoms = mol.atoms, bonds = mol.bonds;
    if (!bonds.some(function (b) { return b.arom; })) return;
    var need = atoms.map(function (a, idx) {
      if (!a.arom) return false;
      var hasExoDouble = mol.adj[idx].some(function (e) { var b = bonds[e.bond]; return !b.arom && b.order === 2; });
      if (hasExoDouble) return false;
      if (a.el === 'C') return true;
      if ((a.el === 'N' || a.el === 'P') && !(a.bracket && a.hExp > 0) && mol.adj[idx].length === 2) return true;
      return false; // o, s, [nH]
    });
    var done = atoms.map(function () { return false; });

    function solve() {
      var k = -1;
      for (var t = 0; t < atoms.length; t++) if (need[t] && !done[t]) { k = t; break; }
      if (k < 0) return true;
      var list = mol.adj[k];
      for (var q = 0; q < list.length; q++) {
        var b = bonds[list[q].bond], nb = list[q].nb;
        if (!b.arom || !need[nb] || done[nb]) continue;
        b.order = 2; done[k] = done[nb] = true;
        if (solve()) return true;
        b.order = 1; done[k] = done[nb] = false;
      }
      return false;
    }
    if (!solve()) throw new Error('Could not kekulize ' + mol.smiles);
  }

  function implicitH(mol) {
    mol.atoms.forEach(function (a, idx) {
      if (a.bracket) { a.h = a.hExp; return; }
      var sum = 0;
      mol.adj[idx].forEach(function (e) { sum += mol.bonds[e.bond].order; });
      var vals = VALENCE[a.el] || [sum];
      var v = vals[vals.length - 1];
      for (var k = 0; k < vals.length; k++) if (vals[k] >= sum) { v = vals[k]; break; }
      a.h = Math.max(0, v - sum);
    });
  }

  /* Smallest ring through each bond (good enough SSSR for small molecules,
     including simple fused systems). */
  function findRings(mol) {
    var rings = [], keys = {};
    mol.bonds.forEach(function (bd, k) {
      var path = shortestPath(mol, bd.a, bd.b, k);
      if (!path) return;
      var key = path.slice().sort(function (x, y) { return x - y; }).join(',');
      if (keys[key]) return;
      keys[key] = true;
      rings.push(path);
    });
    mol.rings = rings;
    mol.atomRings = mol.atoms.map(function () { return []; });
    rings.forEach(function (r, ri) { r.forEach(function (a) { mol.atomRings[a].push(ri); }); });
    mol.bonds.forEach(function (b) {
      b.ring = -1;
      for (var ri = 0; ri < rings.length; ri++) {
        var r = rings[ri], n = r.length;
        for (var k = 0; k < n; k++) {
          var u = r[k], v = r[(k + 1) % n];
          if ((u === b.a && v === b.b) || (u === b.b && v === b.a)) { b.ring = ri; break; }
        }
        if (b.ring >= 0) break;
      }
    });
  }

  function shortestPath(mol, from, to, skipBond) {
    var prev = {}; prev[from] = -1;
    var q = [from];
    while (q.length) {
      var x = q.shift();
      if (x === to) break;
      mol.adj[x].forEach(function (e) {
        if (e.bond === skipBond || prev.hasOwnProperty(e.nb)) return;
        prev[e.nb] = x; q.push(e.nb);
      });
    }
    if (!prev.hasOwnProperty(to)) return null;
    var path = [], cur = to;
    while (cur !== -1) { path.unshift(cur); cur = prev[cur]; }
    return path;
  }

  /* ---------------- 2D layout ---------------- */

  function unit(angDeg) { return { x: Math.cos(angDeg * DEG), y: Math.sin(angDeg * DEG) }; }
  function angleOf(dx, dy) { return Math.atan2(dy, dx) / DEG; }

  function isLinear(mol, idx) {
    var doubles = 0, triple = false;
    mol.adj[idx].forEach(function (e) {
      var o = mol.bonds[e.bond].order;
      if (o === 3) triple = true;
      if (o === 2) doubles++;
    });
    return triple || doubles >= 2;
  }

  /* Longest path (atom count) starting at `from`, not stepping back to `avoid`. */
  function reach(mol, from, avoid) {
    var best = 0, seen = {};
    seen[avoid] = true;
    (function dfs(x, d) {
      seen[x] = true;
      if (d > best) best = d;
      mol.adj[x].forEach(function (e) { if (!seen[e.nb]) dfs(e.nb, d + 1); });
      seen[x] = false;
    })(from, 1);
    return best;
  }

  function bfsDist(mol, from) {
    var dist = {}; dist[from] = 0;
    var q = [from];
    while (q.length) {
      var x = q.shift();
      mol.adj[x].forEach(function (e) { if (!dist.hasOwnProperty(e.nb)) { dist[e.nb] = dist[x] + 1; q.push(e.nb); } });
    }
    return dist;
  }

  function layout(mol) {
    var atoms = mol.atoms, n = atoms.length, L = CFG.bond;
    var pos = new Array(n), st = new Array(n);
    var ringPlaced = mol.rings.map(function () { return false; });
    var queue = [];

    function place(idx, x, y, info) { pos[idx] = { x: x, y: y }; st[idx] = info; queue.push(idx); }

    if (n === 1) { pos[0] = { x: 0, y: 0 }; mol.pos = pos; return; }

    // Start at the end of the longest chain (prefer carbon); ring-only molecules start in a ring.
    var terminals = [];
    for (var t = 0; t < n; t++) if (mol.adj[t].length === 1) terminals.push(t);
    var start = -1;
    if (terminals.length) {
      var bestScore = -1;
      terminals.forEach(function (tt) {
        var d = bfsDist(mol, tt), ecc = 0;
        for (var key in d) if (d[key] > ecc) ecc = d[key];
        var score = ecc * 10 + (atoms[tt].el === 'C' ? 1 : 0);
        if (score > bestScore) { bestScore = score; start = tt; }
      });
    }

    if (start >= 0 && !mol.atomRings[start].length) {
      place(start, 0, 0, { theta: null, lastTurn: -1, parent: -1 });
    } else {
      // pick a ring atom carrying a substituent (if any) so the substituent can leave at -30°
      var r0 = -1;
      for (var a = 0; a < n && r0 < 0; a++) {
        if (!mol.atomRings[a].length) continue;
        if (mol.adj[a].some(function (e) { return !mol.atomRings[e.nb].length || mol.bonds[e.bond].ring < 0; })) r0 = a;
      }
      if (r0 < 0) r0 = mol.rings[0][0];
      var hasExo = mol.adj[r0].some(function (e) { return mol.bonds[e.bond].ring < 0; });
      pos[r0] = { x: 0, y: 0 };
      placeRingSystem(r0, hasExo ? 150 : 90);
    }

    function placeRingSystem(entry, arrive) {
      // first ring: the one containing `entry` (prefer 6-rings)
      var cands = mol.atomRings[entry].slice().sort(function (p, q) {
        return Math.abs(mol.rings[p].length - 6) - Math.abs(mol.rings[q].length - 6);
      });
      var ri = cands[0], ring = mol.rings[ri], m = ring.length;
      var rad = L / (2 * Math.sin(Math.PI / m));
      var u = unit(arrive);
      var cx = pos[entry].x + rad * u.x, cy = pos[entry].y + rad * u.y;
      var k0 = ring.indexOf(entry);
      var a0 = angleOf(pos[entry].x - cx, pos[entry].y - cy);
      var placedRing = [];
      for (var k = 0; k < m; k++) {
        var at = ring[(k0 + k) % m];
        var p = unit(a0 + k * 360 / m);
        if (at !== entry) pos[at] = { x: cx + rad * p.x, y: cy + rad * p.y };
        placedRing.push(at);
      }
      ringPlaced[ri] = true;
      var centers = {}; centers[ri] = { x: cx, y: cy };

      // fused rings sharing an edge with an already-placed ring
      var progress = true;
      while (progress) {
        progress = false;
        for (var rj = 0; rj < mol.rings.length; rj++) {
          if (ringPlaced[rj]) continue;
          var rr = mol.rings[rj], mm = rr.length;
          var shared = rr.filter(function (x) { return pos[x] && placedRing.indexOf(x) >= 0; });
          if (shared.length !== 2) continue;
          var su = shared[0], sv = shared[1];
          if (!bondBetween(mol, su, sv)) continue;
          // which placed ring owns this edge?
          var owner = null;
          for (var key in centers) {
            var ring2 = mol.rings[key];
            if (ring2.indexOf(su) >= 0 && ring2.indexOf(sv) >= 0) { owner = centers[key]; break; }
          }
          var mx = (pos[su].x + pos[sv].x) / 2, my = (pos[su].y + pos[sv].y) / 2;
          var ex = pos[sv].x - pos[su].x, ey = pos[sv].y - pos[su].y, el = Math.sqrt(ex * ex + ey * ey);
          var nx = -ey / el, ny = ex / el;
          if (owner && ((owner.x - mx) * nx + (owner.y - my) * ny) > 0) { nx = -nx; ny = -ny; }
          var apo = L / (2 * Math.tan(Math.PI / mm)), rad2 = L / (2 * Math.sin(Math.PI / mm));
          var ncx = mx + apo * nx, ncy = my + apo * ny;
          // walk the new ring from su, in the direction that goes su -> sv first, then reverse to continue past su
          var ku = rr.indexOf(su), kv = rr.indexOf(sv);
          var stepDir = ((ku + 1) % mm === kv) ? 1 : -1;
          var au = angleOf(pos[su].x - ncx, pos[su].y - ncy);
          var av = angleOf(pos[sv].x - ncx, pos[sv].y - ncy);
          var dAng = ((av - au + 540) % 360) - 180; // signed step su -> sv
          for (var s2 = 2; s2 < mm; s2++) {
            var at2 = rr[((ku + stepDir * s2) % mm + mm) % mm];
            var pp = unit(au + dAng * s2);
            pos[at2] = { x: ncx + rad2 * pp.x, y: ncy + rad2 * pp.y };
            placedRing.push(at2);
          }
          ringPlaced[rj] = true;
          centers[rj] = { x: ncx, y: ncy };
          progress = true;
        }
      }
      placedRing.forEach(function (at3) {
        // centre used for exocyclic direction = mean of the centres of the rings this atom sits in
        var sx = 0, sy = 0, c = 0;
        mol.atomRings[at3].forEach(function (q) { if (centers[q]) { sx += centers[q].x; sy += centers[q].y; c++; } });
        var info = { ring: true, cx: sx / c, cy: sy / c };
        st[at3] = info;
        queue.push(at3);
      });
    }

    function placeChild(parent, child, ang, lastTurn) {
      var u = unit(ang);
      var x = pos[parent].x + L * u.x, y = pos[parent].y + L * u.y;
      if (mol.atomRings[child].length && !ringPlaced[mol.atomRings[child][0]]) {
        pos[child] = { x: x, y: y };
        placeRingSystem(child, ang);
      } else {
        place(child, x, y, { theta: ang, lastTurn: lastTurn, parent: parent });
      }
    }

    function priority(A, kid) {
      return reach(mol, kid, A) * 10 + (atoms[kid].el === 'C' ? 1 : 0);
    }

    while (queue.length) {
      var A = queue.shift(), info = st[A];
      var kids = mol.adj[A].map(function (e) { return e.nb; }).filter(function (k) { return !pos[k]; });
      if (!kids.length) continue;
      kids.sort(function (p, q) { return priority(A, q) - priority(A, p); });

      if (info.ring) {
        var radial = angleOf(pos[A].x - info.cx, pos[A].y - info.cy);
        if (kids.length === 1) placeChild(A, kids[0], radial, -1);
        else {
          var spread = kids.length === 2 ? 30 : 45;
          kids.forEach(function (k, ix) { placeChild(A, k, radial - spread + ix * (2 * spread / (kids.length - 1)), -1); });
        }
        continue;
      }

      var th = info.theta, s = -info.lastTurn || 1;
      if (info.parent < 0) {
        kids.forEach(function (k, ix) { placeChild(A, k, -30 + ix * 360 / kids.length, -1); });
      } else if (isLinear(mol, A)) {
        kids.forEach(function (k) { placeChild(A, k, th, info.lastTurn); });
      } else if (kids.length === 1) {
        placeChild(A, kids[0], th + 60 * s, s);
      } else if (kids.length === 2) {
        placeChild(A, kids[0], th + 60 * s, s);
        placeChild(A, kids[1], th - 60 * s, -s);
      } else {
        placeChild(A, kids[0], th + 60 * s, s);
        placeChild(A, kids[1], th - 20 * s, -s);
        placeChild(A, kids[2], th - 100 * s, s);
      }
    }

    mol.pos = pos;
    enforceStereo(mol);
    normaliseOrientation(mol);
  }

  /* E/Z: compare the SMILES-specified relationship with the drawn one and,
     if wrong, reflect one side of the double bond across its axis. */
  function enforceStereo(mol) {
    var pos = mol.pos;
    mol.bonds.forEach(function (bd, k) {
      if (bd.order !== 2 || bd.ring >= 0) return;
      var ra = stereoRef(mol, bd.a, k), rb = stereoRef(mol, bd.b, k);
      if (!ra || !rb) return;
      var wantCis = ra.side === rb.side;
      var A = pos[bd.a], B = pos[bd.b];
      var dx = B.x - A.x, dy = B.y - A.y;
      var ca = dx * (pos[ra.nb].y - A.y) - dy * (pos[ra.nb].x - A.x);
      var cb = dx * (pos[rb.nb].y - B.y) - dy * (pos[rb.nb].x - B.x);
      var isCis = (ca > 0) === (cb > 0);
      if (isCis === wantCis) return;
      // reflect everything on B's side across the A–B line
      var side = {}, q = [bd.b]; side[bd.b] = true;
      while (q.length) {
        var x = q.shift();
        mol.adj[x].forEach(function (e) {
          if (e.bond === k || side[e.nb]) return;
          side[e.nb] = true; q.push(e.nb);
        });
      }
      if (side[bd.a]) return; // in a ring system — leave it
      var len2 = dx * dx + dy * dy;
      Object.keys(side).forEach(function (key) {
        var P = pos[key];
        var t = ((P.x - A.x) * dx + (P.y - A.y) * dy) / len2;
        var fx = A.x + t * dx, fy = A.y + t * dy;
        pos[key] = { x: 2 * fx - P.x, y: 2 * fy - P.y };
      });
    });
  }

  function stereoRef(mol, atom, skipBond) {
    var list = mol.adj[atom];
    for (var q = 0; q < list.length; q++) {
      if (list[q].bond === skipBond) continue;
      var b = mol.bonds[list[q].bond];
      if (!b.dir) continue;
      return { nb: list[q].nb, side: (b.a === atom) ? b.dir : -b.dir };
    }
    return null;
  }

  /* Rotate in 30° steps (keeps the hexagonal grid) to the flattest orientation.
     Skipped for molecules with linear sp centres: a straight C–C≡C–C run
     already sits on the textbook 30° diagonal, and turning it horizontal
     makes the rest of the chain kink at 60° instead of zigzagging. */
  function normaliseOrientation(mol) {
    for (var i = 0; i < mol.atoms.length; i++) if (isLinear(mol, i)) return;
    var pos = mol.pos, best = null;
    for (var k = 0; k < 12; k++) {
      var ang = k * 30 * DEG, c = Math.cos(ang), s = Math.sin(ang);
      var minY = Infinity, maxY = -Infinity;
      pos.forEach(function (p) {
        var y = p.x * s + p.y * c;
        if (y < minY) minY = y; if (y > maxY) maxY = y;
      });
      var score = (maxY - minY) + (k === 0 ? -0.5 : 0); // slight preference to leave as drawn
      if (!best || score < best.score - 1e-6) best = { score: score, c: c, s: s };
    }
    mol.pos = pos.map(function (p) { return { x: p.x * best.c - p.y * best.s, y: p.x * best.s + p.y * best.c }; });
  }

  /* ---------------- Rendering ---------------- */

  var CW = { H: 0.72, C: 0.72, N: 0.72, O: 0.78, S: 0.67, P: 0.67, B: 0.67, F: 0.61, I: 0.28, l: 0.22, r: 0.33, '+': 0.58, '\u2212': 0.58 };
  function textWidth(str, size) {
    var w = 0;
    for (var i = 0; i < str.length; i++) w += (CW[str[i]] || (/[0-9]/.test(str[i]) ? 0.56 : 0.6));
    return w * size;
  }

  function f(v) { return (Math.round(v * 10) / 10).toString(); }

  function needsLabel(mol, idx) {
    var a = mol.atoms[idx];
    return a.el !== 'C' || a.charge !== 0 || mol.adj[idx].length === 0;
  }

  function labelGeometry(mol, idx) {
    var a = mol.atoms[idx], P = mol.pos[idx], F = CFG.font;
    var symW = textWidth(a.el, F);
    var hStr = a.h ? 'H' + (a.h > 1 ? a.h : '') : '';
    var hW = a.h ? textWidth('H', F) + (a.h > 1 ? textWidth(String(a.h), F * 0.75) : 0) : 0;
    var sumDx = 0;
    mol.adj[idx].forEach(function (e) { sumDx += mol.pos[e.nb].x - P.x; });
    var hLeft = sumDx > 1; // bonds come from the right -> write HO, H2N
    return {
      symW: symW, hW: hW, hLeft: hLeft, h: a.h,
      hw: symW / 2 + 1.5, hh: F * 0.46,
      left: P.x - symW / 2 - (hLeft ? hW : 0),
      right: P.x + symW / 2 + (hLeft ? 0 : hW) + (a.charge ? F * 0.5 : 0)
    };
  }

  function clipDist(g, ux, uy) {
    if (!g) return 0;
    var ax = Math.abs(ux), ay = Math.abs(uy);
    var dx = ax > 1e-6 ? g.hw / ax : Infinity;
    var dy = ay > 1e-6 ? g.hh / ay : Infinity;
    return Math.min(dx, dy) + 1;
  }

  function line(x1, y1, x2, y2) {
    return '<line x1="' + f(x1) + '" y1="' + f(y1) + '" x2="' + f(x2) + '" y2="' + f(y2) + '"/>';
  }

  function drawBonds(mol, labels) {
    var pos = mol.pos, L = CFG.bond, out = '';
    mol.bonds.forEach(function (bd) {
      var A = pos[bd.a], B = pos[bd.b];
      var dx = B.x - A.x, dy = B.y - A.y, len = Math.sqrt(dx * dx + dy * dy);
      var ux = dx / len, uy = dy / len, nx = -uy, ny = ux;
      var ca = clipDist(labels[bd.a], ux, uy), cb = clipDist(labels[bd.b], ux, uy);
      var x1 = A.x + ux * ca, y1 = A.y + uy * ca, x2 = B.x - ux * cb, y2 = B.y - uy * cb;
      var g = CFG.gap * L;

      if (bd.order === 1) { out += line(x1, y1, x2, y2); return; }
      if (bd.order === 3) {
        out += line(x1, y1, x2, y2);
        out += line(x1 + nx * g, y1 + ny * g, x2 + nx * g, y2 + ny * g);
        out += line(x1 - nx * g, y1 - ny * g, x2 - nx * g, y2 - ny * g);
        return;
      }
      // double bond: decide which side the second line goes
      var mx = (A.x + B.x) / 2, my = (A.y + B.y) / 2, side = 0, centred = false;
      if (bd.ring >= 0) {
        var r = mol.rings[bd.ring], cx = 0, cy = 0;
        r.forEach(function (at) { cx += pos[at].x; cy += pos[at].y; });
        cx /= r.length; cy /= r.length;
        side = ((cx - mx) * nx + (cy - my) * ny) > 0 ? 1 : -1;
      } else if (mol.adj[bd.a].length === 1 || mol.adj[bd.b].length === 1) {
        centred = true;
      } else {
        var score = 0;
        [bd.a, bd.b].forEach(function (end) {
          mol.adj[end].forEach(function (e) {
            if (e.nb === bd.a || e.nb === bd.b) return;
            var d = (pos[e.nb].x - mx) * nx + (pos[e.nb].y - my) * ny;
            score += d > 0.5 ? 1 : d < -0.5 ? -1 : 0;
          });
        });
        if (score === 0) centred = true; else side = score > 0 ? 1 : -1;
      }
      if (centred) {
        var h = g / 2;
        out += line(x1 + nx * h, y1 + ny * h, x2 + nx * h, y2 + ny * h);
        out += line(x1 - nx * h, y1 - ny * h, x2 - nx * h, y2 - ny * h);
      } else {
        var ta = labels[bd.a] ? 0 : CFG.trim * L, tb = labels[bd.b] ? 0 : CFG.trim * L;
        out += line(x1, y1, x2, y2);
        out += line(x1 + ux * ta + nx * g * side, y1 + uy * ta + ny * g * side,
                    x2 - ux * tb + nx * g * side, y2 - uy * tb + ny * g * side);
      }
    });
    return out;
  }

  function drawLabels(mol, labels) {
    var F = CFG.font, sub = F * 0.75, out = '';
    labels.forEach(function (g, idx) {
      if (!g) return;
      var a = mol.atoms[idx], P = mol.pos[idx];
      var base = P.y + F * 0.36;
      out += '<text x="' + f(P.x) + '" y="' + f(base) + '" text-anchor="middle">' + a.el + '</text>';
      if (g.h) {
        var hTxt = 'H' + (g.h > 1 ? '<tspan dy="' + f(F * 0.28) + '" font-size="' + f(sub) + '">' + g.h + '</tspan>' : '');
        if (g.hLeft) out += '<text x="' + f(P.x - g.symW / 2) + '" y="' + f(base) + '" text-anchor="end">' + hTxt + '</text>';
        else out += '<text x="' + f(P.x + g.symW / 2) + '" y="' + f(base) + '" text-anchor="start">' + hTxt + '</text>';
      }
      if (a.charge) {
        var q = Math.abs(a.charge);
        var chTxt = (q > 1 ? q : '') + (a.charge > 0 ? '+' : '\u2212');
        var cx = (g.hLeft ? P.x + g.symW / 2 : g.right - F * 0.5) + 0.5;
        out += '<text x="' + f(cx) + '" y="' + f(P.y - F * 0.2) + '" font-size="' + f(sub) + '" text-anchor="start">' + chTxt + '</text>';
      }
    });
    return out;
  }

  function prepare(smiles) {
    var mol = parse(smiles);
    layout(mol);
    var labels = mol.atoms.map(function (a, idx) { return needsLabel(mol, idx) ? labelGeometry(mol, idx) : null; });
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    mol.pos.forEach(function (p, idx) {
      var g = labels[idx];
      var l = g ? g.left : p.x, r = g ? g.right : p.x;
      var t = g ? p.y - g.hh - (mol.atoms[idx].charge ? 3 : 0) : p.y;
      var b = g ? p.y + g.hh + (g.h > 1 ? 3 : 0) : p.y;
      if (l < minX) minX = l; if (r > maxX) maxX = r;
      if (t < minY) minY = t; if (b > maxY) maxY = b;
    });
    var M = CFG.margin + CFG.lineWidth;
    return {
      mol: mol, labels: labels,
      box: { x: minX - M, y: minY - M, w: (maxX - minX) + 2 * M, h: (maxY - minY) + 2 * M }
    };
  }

  function toSVG(p, opts) {
    opts = opts || {};
    var scale = opts.scale || 1;
    var box = p.box, w = box.w, h = box.h, y0 = box.y, x0 = box.x;
    if (opts.minHeight && opts.minHeight > h) { y0 -= (opts.minHeight - h) / 2; h = opts.minHeight; }
    if (opts.minWidth && opts.minWidth > w) { x0 -= (opts.minWidth - w) / 2; w = opts.minWidth; }
    return '<svg xmlns="http://www.w3.org/2000/svg" class="mol" role="img"' +
      (opts.title ? ' aria-label="' + opts.title + '"' : '') +
      ' viewBox="' + f(x0) + ' ' + f(y0) + ' ' + f(w) + ' ' + f(h) + '"' +
      ' width="' + Math.round(w * scale) + '" height="' + Math.round(h * scale) + '">' +
      '<g fill="none" stroke="currentColor" stroke-width="' + CFG.lineWidth + '" stroke-linecap="round">' +
      drawBonds(p.mol, p.labels) + '</g>' +
      '<g fill="currentColor" stroke="none" font-family="' + CFG.fontFamily + '" font-size="' + CFG.font + '">' +
      drawLabels(p.mol, p.labels) + '</g></svg>';
  }

  function svg(smiles, opts) { return toSVG(prepare(smiles), opts); }

  /* Several structures drawn at the same scale and padded to a common height,
     so captions above them line up. */
  function row(smilesList, opts) {
    var ps = smilesList.map(prepare);
    var hMax = Math.max.apply(null, ps.map(function (p) { return p.box.h; }));
    return ps.map(function (p) {
      var o = {}; for (var k in (opts || {})) o[k] = opts[k];
      o.minHeight = hMax;
      return toSVG(p, o);
    });
  }

  return { CFG: CFG, parse: parse, svg: svg, row: row, pair: function (a, b, o) { return row([a, b], o); } };
})();
