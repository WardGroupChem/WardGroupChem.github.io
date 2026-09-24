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

/* ============================================================
   STRUCTURES — drawn by chem.js from SMILES (see that file).
   Every molecule uses the same bond length, angles and label rules.
   ============================================================ */

function renderReaction(startSmiles, targetSmiles) {
  var svgs;
  try {
    svgs = Chem.pair(startSmiles, targetSmiles);
  } catch (e) {
    // Fallback so a drawing problem can never freeze the game/quiz
    if (window.console) console.error('Structure drawing failed:', e);
    svgs = ['<code>' + startSmiles + '</code>', '<code>' + targetSmiles + '</code>'];
  }
  var html = '<div class="reaction-row">';
  html += '<div class="reaction-mol"><span class="mol-cap">Start</span>' + svgs[0] + '</div>';
  html += '<div class="reaction-arrow"><span class="arrow-q">?</span><span class="arrow-line"><span class="arrow-shaft"></span><span class="arrow-head"></span></span></div>';
  html += '<div class="reaction-mol"><span class="mol-cap">Target</span>' + svgs[1] + '</div>';
  html += '</div>';
  return html;
}

/* ============================================================
   PUZZLE DATA — alkene & alkyne synthesis, 2nd year undergrad
   ============================================================ */

var PUZZLES = [
  {
    title: 'But-2-ene, from an alkyl halide',
    start: 'CC(Br)CC',
    target: 'C/C=C/C',
    options: [
      {
        correct: true,
        reagents: 'NaOEt, EtOH, Δ',
        explain: 'Correct — a small, strong base like ethoxide favours clean E2 elimination to the more substituted (Zaitsev) alkene, but-2-ene, as the major product.'
      },
      {
        correct: false,
        reagents: 'KOtBu, Δ',
        explain: 'A bulky base like potassium tert-butoxide can\'t easily reach the more hindered internal proton, so it removes the more accessible terminal one instead — favouring the less substituted (Hofmann) alkene, but-1-ene, not the product shown.'
      },
      {
        correct: false,
        reagents: 'NaCN, DMSO',
        explain: 'Cyanide is a good nucleophile in a polar aprotic solvent — this favours SN2 substitution to the nitrile rather than elimination to an alkene.'
      },
      {
        correct: false,
        reagents: 'Br₂, CCl₄',
        explain: 'Br₂ adds across a C=C double bond — it\'s a test/trap for an alkene, not a reagent that installs one from an alkyl halide.'
      }
    ]
  },
  {
    title: 'But-2-ene, from an alcohol',
    start: 'CC(O)CC',
    target: 'C/C=C/C',
    options: [
      {
        correct: true,
        reagents: 'H₂SO₄ (conc.), Δ',
        explain: 'Correct — acid protonates the OH into a good leaving group (water); E1 loss of water via the secondary carbocation, then loss of the more substituted β-H, gives the more stable (Zaitsev) alkene, but-2-ene, as the major product.'
      },
      {
        correct: false,
        reagents: 'PCC, CH₂Cl₂',
        explain: 'PCC oxidises this secondary alcohol to the ketone, butan-2-one — a useful reaction, but not a dehydration, so no alkene forms.'
      },
      {
        correct: false,
        reagents: 'SOCl₂, pyridine',
        explain: 'This swaps the OH for Cl, giving the alkyl chloride — no elimination happens without a separate base-mediated step afterwards.'
      },
      {
        correct: false,
        reagents: 'NaOH (aq)',
        explain: 'Hydroxide can\'t protonate the OH into a leaving group — without acid activation, this alcohol won\'t eliminate under basic conditions.'
      }
    ]
  },
  {
    title: 'Stilbene, by Wittig olefination',
    start: 'O=Cc1ccccc1',
    target: 'c1ccc(cc1)C=Cc1ccccc1',
    options: [
      {
        correct: true,
        reagents: 'PPh₃, PhCH₂Br, then n-BuLi',
        explain: 'Correct — PPh₃ and benzyl bromide form a phosphonium salt; deprotonation with base generates the benzylidene ylide, which reacts with the aldehyde shown to form the new C=C bond directly (a Wittig reaction). This semi-stabilised ylide typically gives a mixture of (E)- and (Z)-stilbene.'
      },
      {
        correct: false,
        reagents: 'PhMgBr, then H₃O⁺',
        explain: 'A Grignard reagent adds to the aldehyde to give a secondary alcohol (1,2-diphenylethan-1-ol) after work-up — an alkene never forms, since no ylide is involved.'
      },
      {
        correct: false,
        reagents: 'H₂NNH₂, KOH, Δ',
        explain: 'Wolff–Kishner conditions simply reduce the aldehyde to a methyl group — no new carbon–carbon bond is formed to the second ring.'
      },
      {
        correct: false,
        reagents: 'NaBH₄',
        explain: 'Sodium borohydride reduces the aldehyde straight to benzyl alcohol — again, no new C–C bond is formed.'
      }
    ]
  },
  {
    title: '(Z)-Pent-2-ene, from an alkyne',
    start: 'CC#CCC',
    target: 'C/C=C\\CC',
    options: [
      {
        correct: true,
        reagents: 'H₂, Lindlar catalyst',
        explain: 'Correct — Lindlar\'s poisoned Pd catalyst delivers H₂ syn across the triple bond in a single step and stops at the alkene, giving the cis (Z) product.'
      },
      {
        correct: false,
        reagents: 'Na, NH₃ (l)',
        explain: 'Dissolving-metal reduction does reduce the alkyne to an alkene, but via anti addition — it gives the trans (E) isomer, not the cis (Z) alkene shown.'
      },
      {
        correct: false,
        reagents: 'H₂, Pd/C (excess)',
        explain: 'An unpoisoned Pd/C catalyst doesn\'t stop at the alkene stage — it keeps reducing all the way to the fully saturated alkane, pentane.'
      },
      {
        correct: false,
        reagents: 'HBr (1 eq)',
        explain: 'HBr adds H and Br across the triple bond (Markovnikov) to give a bromoalkene — it doesn\'t simply hydrogenate to the alkene shown.'
      }
    ]
  },
  {
    title: '(E)-Pent-2-ene, from an alkyne',
    start: 'CC#CCC',
    target: 'C/C=C/CC',
    options: [
      {
        correct: true,
        reagents: 'Na, NH₃ (l)',
        explain: 'Correct — dissolving-metal reduction proceeds through a trans-vinyl radical intermediate, delivering the two new C–H bonds anti to each other and giving the trans (E) alkene.'
      },
      {
        correct: false,
        reagents: 'H₂, Lindlar catalyst',
        explain: 'Lindlar\'s catalyst delivers syn hydrogenation, giving the cis (Z) alkene — the opposite geometry to the product shown.'
      },
      {
        correct: false,
        reagents: 'H₂, Pd/C (excess)',
        explain: 'An unpoisoned Pd/C catalyst over-reduces straight through to the alkane, pentane, rather than stopping at the alkene.'
      },
      {
        correct: false,
        reagents: 'HBr (2 eq)',
        explain: 'Two equivalents of HBr add twice across the triple bond, giving a geminal dihalide rather than stopping at the alkene stage.'
      }
    ]
  },
  {
    title: 'But-1-yne, from a vicinal dihalide',
    start: 'BrCC(Br)CC',
    target: 'C#CCC',
    options: [
      {
        correct: true,
        reagents: 'LDA (excess), THF, Δ; then H₂O',
        explain: 'Correct — excess strong, hindered, non-nucleophilic base drives two successive E2 eliminations (via the bromoalkene) to the alkyne, and also deprotonates it to the acetylide; aqueous work-up then protonates it back to but-1-yne.'
      },
      {
        correct: false,
        reagents: 'KOH, EtOH, Δ',
        explain: 'Hydroxide is strong enough for the first elimination but isn\'t up to the second, far less acidic vinylic proton — the reaction typically stalls at the bromoalkene stage.'
      },
      {
        correct: false,
        reagents: 'H₂, Pd/C',
        explain: 'This is a hydrogenation catalyst — there\'s nothing here to remove from an alkyl halide, and it won\'t form the alkyne.'
      },
      {
        correct: false,
        reagents: 'Mg, Et₂O',
        explain: 'These are Grignard-forming conditions (making an organometallic from the C–Br bond) — a completely different reaction type, not an elimination.'
      }
    ]
  },
  {
    title: 'Propyne, from a geminal dihalide',
    start: 'CC(Br)(Br)C',
    target: 'CC#C',
    options: [
      {
        correct: true,
        reagents: 'LDA (excess), THF, Δ; then H₂O',
        explain: 'Correct — excess strong, hindered base removes both C–Br bonds from the same carbon in two successive E2 steps (via the bromoalkene), and aqueous work-up protonates the resulting acetylide to give propyne.'
      },
      {
        correct: false,
        reagents: 'NaOH (aq), Δ',
        explain: 'Hydroxide can manage the first elimination but isn\'t strong enough for the second, far less acidic vinylic proton — the reaction typically stalls at 2-bromopropene.'
      },
      {
        correct: false,
        reagents: 'Zn, AcOH',
        explain: 'Zinc in acetic acid reductively removes halogens (e.g. debrominating a vicinal dihalide back to an alkene) — it isn\'t a base, and won\'t drive elimination on to the alkyne.'
      },
      {
        correct: false,
        reagents: 'AgNO₃, EtOH',
        explain: 'Silver nitrate promotes ionisation of the C–Br bonds (SN1-type solvolysis) rather than the E2 double elimination needed to reach the alkyne.'
      }
    ]
  },
  {
    title: 'Hept-3-yne, by acetylide alkylation',
    start: 'C#CCCC',
    target: 'CCC#CCCC',
    options: [
      {
        correct: true,
        reagents: 'LDA, THF; then CH₃CH₂Br',
        explain: 'Correct — LDA deprotonates the terminal alkyne to give the acetylide, a strong carbanion nucleophile, which performs SN2 substitution on the primary ethyl bromide to extend the chain.'
      },
      {
        correct: false,
        reagents: 'LDA, THF; then (CH₃)₃CBr',
        explain: 'The acetylide is a strong base — with a tertiary halide like tert-butyl bromide it promotes E2 elimination instead of SN2 substitution, so the chain isn\'t extended as drawn.'
      },
      {
        correct: false,
        reagents: 'CH₃CH₂Br (no base first)',
        explain: 'Without first forming the acetylide, there\'s no nucleophile present to attack the alkyl bromide — nothing happens.'
      },
      {
        correct: false,
        reagents: 'NaOH; then CH₃CH₂Br',
        explain: 'Hydroxide (conjugate acid pKₐ ~15.7) is nowhere near basic enough to deprotonate a terminal alkyne (pKₐ ~25) — no acetylide forms, so no alkylation can occur.'
      }
    ]
  }
];

/* ============================================================
   INVADER SPRITES — classic blocky 8-bit alien silhouettes,
   recoloured into the project palette (chemistry stays in the
   bullets/capsule; the aliens stay genuinely "space invader").
   Two 5x6 pixel-grid variants, alternated by row.
   ============================================================ */
var INV_PIXEL = 2;
var SPRITES = [
  ['.XXX.',
   'XXXXX',
   'XX.XX',
   'XXXXX',
   '.X.X.',
   'X...X'],
  ['X...X',
   '.X.X.',
   'XXXXX',
   'XX.XX',
   'XXXXX',
   '.X.X.']
];
var SPRITE_W = 5 * INV_PIXEL, SPRITE_H = 6 * INV_PIXEL;

/* ============================================================
   CONSTANTS
   ============================================================ */
var CANVAS_W = 384, CANVAS_H = 480;
var INV_COLS = 8, INV_ROWS = 5;
var INV_PITCH_X = 40, INV_PITCH_Y = 30;
var INV_ORIGIN_X = 32, INV_ORIGIN_Y = 40;
var INV_BASE_SPEED = 0.2;
var INV_WAVE_STEP = 0.06;
var INV_DROP_STEP = 10;
var PLAYER_SPEED = 3.2;
var PLAYER_BULLET_SPEED = 4.2; // was 6 — 30% slower so the bond labels are readable
var FIRE_COOLDOWN_FRAMES = 14;
var CAPSULE_INTERVAL_FRAMES = 1200; // ~20s at 60fps
var SHIELD_INTERVAL_FRAMES = 600; // ~10s at 60fps
var SHIELD_FALL_SPEED = 0.9;
var SHIELD_COLOR = '#4FD8E8';
var TIER_LABELS = ['Alkane —', 'Alkene =', 'Alkyne \u2261'];
var TIER_COLORS = ['#8A8FA8', '#6C5CE7', '#00BFA6'];

var ALIEN_FIRE_CHECK_FRAMES = 90;
var ALIEN_FIRE_CHANCE = 0.22;
var ALIEN_BULLET_SPEED = 2.2;
var BULLET_LABELS = ['Me', 'CH\u2082', 'CH'];

var INCOMING_LEAD_FRAMES = 60; // ~1s warning before a capsule/shield drop
var WAVE_BANNER_FRAMES = 90;
var LIFE_BANNER_FRAMES = 60;
var SHAKE_FRAMES = 14;
var SHAKE_MAGNITUDE = 6;

var DIFFICULTY_PRESETS = {
  easy:   { lives: 7, fireChance: 0.13, waveStep: 0.045 },
  normal: { lives: 5, fireChance: 0.22, waveStep: 0.06 },
  hard:   { lives: 3, fireChance: 0.32, waveStep: 0.08 }
};
var difficulty = 'normal';

var BUNKER_TEMPLATE = [
  '...XXXXX...',
  '..XXXXXXX..',
  '.XXXXXXXXX.',
  'XXXXXXXXXXX',
  'XXXXXXXXXXX',
  '.XXXXXXXXX.',
  '..XXXXXXX..',
  '...XXXXX...'
];
var BUNKER_PIXEL = 4;
var BUNKER_COLS = 11, BUNKER_ROWS = 8;
var BUNKER_COUNT = 4;

/* ============================================================
   STATE
   ============================================================ */
var game = {
  lives: 5,
  liveScore: 0,
  solved: [false, false, false, false, false, false, false, false],
  picks: [],
  askedCount: 0,
  bulletTier: 0,
  wave: 1,
  quizPaused: false,
  manualPaused: false,
  running: false,
  hasShield: false,
  banner: null
};

var canvas, ctx, player, bullets, invaders, groupX, groupY, groupDir;
var capsule, questionOrder, activeQuestionIdx;
var capsuleState = 'idle', capsuleTimer = 0;
var shieldBonus;
var shieldState = 'idle', shieldTimer = 0;
var incomingCue = null;
var shakeFramesLeft = 0;
var particles = [];
var fireCooldown = 0;
var keys = { left: false, right: false };
var bunkers, invaderBullets, alienFireTimer, frameCount;
var STAR_FIELD = [], BG_HEXES = [];

function $(id) { return document.getElementById(id); }

/* ============================================================
   SFX — tiny synthesized sound effects (oscillators + a noise
   burst), so the package stays self-contained with no audio
   files to bundle or fail to load inside an LMS iframe. Muted by
   default; the mute button in the HUD toggles this.
   ============================================================ */
var SFX = (function () {
  var ctx = null;
  var muted = true;

  function ensureCtx() {
    if (!ctx) {
      var AC = window.AudioContext || window.webkitAudioContext;
      if (AC) { try { ctx = new AC(); } catch (e) { ctx = null; } }
    }
    if (ctx && ctx.state === 'suspended') { try { ctx.resume(); } catch (e2) {} }
    return ctx;
  }

  function tone(freq, dur, type, opts) {
    if (muted) return;
    var c = ensureCtx();
    if (!c) return;
    opts = opts || {};
    var t0 = c.currentTime;
    var osc = c.createOscillator();
    var gain = c.createGain();
    osc.type = type || 'square';
    osc.frequency.setValueAtTime(freq, t0);
    if (opts.slideTo) osc.frequency.exponentialRampToValueAtTime(Math.max(opts.slideTo, 1), t0 + dur);
    var peak = opts.volume != null ? opts.volume : 0.16;
    gain.gain.setValueAtTime(0.0001, t0);
    gain.gain.exponentialRampToValueAtTime(peak, t0 + 0.01);
    gain.gain.exponentialRampToValueAtTime(0.0001, t0 + dur);
    osc.connect(gain); gain.connect(c.destination);
    osc.start(t0);
    osc.stop(t0 + dur + 0.02);
  }

  function noiseBurst(dur, opts) {
    if (muted) return;
    var c = ensureCtx();
    if (!c) return;
    opts = opts || {};
    var t0 = c.currentTime;
    var bufferSize = Math.max(1, Math.floor(c.sampleRate * dur));
    var buffer = c.createBuffer(1, bufferSize, c.sampleRate);
    var data = buffer.getChannelData(0);
    for (var i = 0; i < bufferSize; i++) { data[i] = (Math.random() * 2 - 1) * (1 - i / bufferSize); }
    var src = c.createBufferSource();
    src.buffer = buffer;
    var filter = c.createBiquadFilter();
    filter.type = 'lowpass';
    filter.frequency.setValueAtTime(opts.filterFreq || 1200, t0);
    var gain = c.createGain();
    gain.gain.setValueAtTime(opts.volume || 0.22, t0);
    gain.gain.exponentialRampToValueAtTime(0.0001, t0 + dur);
    src.connect(filter); filter.connect(gain); gain.connect(c.destination);
    src.start(t0);
  }

  return {
    setMuted: function (m) { muted = !!m; },
    isMuted: function () { return muted; },
    unlock: function () { ensureCtx(); },
    fire: function () { tone(620, 0.07, 'square', { slideTo: 340, volume: 0.09 }); },
    pop: function () { noiseBurst(0.14, { filterFreq: 1600, volume: 0.2 }); },
    shield: function () {
      tone(320, 0.12, 'sine', { slideTo: 680, volume: 0.14 });
      setTimeout(function () { tone(680, 0.1, 'sine', { slideTo: 900, volume: 0.1 }); }, 70);
    },
    correct: function () {
      tone(520, 0.09, 'triangle', { volume: 0.15 });
      setTimeout(function () { tone(780, 0.16, 'triangle', { volume: 0.16 }); }, 90);
    },
    wrong: function () { tone(220, 0.22, 'sawtooth', { slideTo: 100, volume: 0.16 }); },
    lifeLost: function () { tone(200, 0.25, 'sawtooth', { slideTo: 55, volume: 0.18 }); },
    waveClear: function () {
      tone(440, 0.1, 'triangle', { volume: 0.13 });
      setTimeout(function () { tone(660, 0.16, 'triangle', { volume: 0.14 }); }, 110);
    }
  };
})();

function shuffle(arr) {
  var a = arr.slice();
  for (var i = a.length - 1; i > 0; i--) {
    var j = Math.floor(Math.random() * (i + 1));
    var t = a[i]; a[i] = a[j]; a[j] = t;
  }
  return a;
}

/* ============================================================
   INVADERS
   ============================================================ */
function initInvaders() {
  invaders = [];
  // Randomise which of the two sprite shapes each row uses, freshly each
  // wave (including the first), so playthroughs don't all look identical.
  var rowVariant = [];
  for (var r = 0; r < INV_ROWS; r++) rowVariant[r] = Math.random() < 0.5 ? 0 : 1;
  for (var r2 = 0; r2 < INV_ROWS; r2++) {
    for (var c = 0; c < INV_COLS; c++) {
      invaders.push({ row: r2, col: c, alive: true, variant: rowVariant[r2] });
    }
  }
  groupX = 0; groupY = 0; groupDir = 1;
}

function invaderPos(inv) {
  return {
    x: Math.round(INV_ORIGIN_X + inv.col * INV_PITCH_X + groupX),
    y: Math.round(INV_ORIGIN_Y + inv.row * INV_PITCH_Y + groupY)
  };
}

function resetWaveAfterHit() {
  invaders.forEach(function (inv) { inv.alive = true; });
  groupX = 0; groupY = 0; groupDir = 1;
  bullets = [];
  invaderBullets = [];
}

function updateInvaders() {
  var alive = invaders.filter(function (i) { return i.alive; });
  if (alive.length === 0) {
    if (!game.banner) {
      var clearedWave = game.wave;
      SFX.waveClear();
      startBanner('Wave ' + clearedWave + ' cleared', 'Wave ' + (clearedWave + 1) + ' incoming…', WAVE_BANNER_FRAMES, 'wave', function () {
        game.wave = clearedWave + 1;
        initInvaders();
      });
    }
    return;
  }

  var frac = alive.length / (INV_ROWS * INV_COLS);
  var speed = (INV_BASE_SPEED + (game.wave - 1) * INV_WAVE_STEP) * (1 + (1 - frac) * 1.6);
  groupX += speed * groupDir;

  var minX = Infinity, maxX = -Infinity, maxY = -Infinity;
  alive.forEach(function (inv) {
    var p = invaderPos(inv);
    if (p.x < minX) minX = p.x;
    if (p.x > maxX) maxX = p.x;
    if (p.y > maxY) maxY = p.y;
  });

  if (maxX + SPRITE_W > CANVAS_W - 10 || minX < 10) {
    groupDir *= -1;
    groupY += INV_DROP_STEP;
  }

  if (maxY + SPRITE_H > player.y - 14) {
    loseLife();
    if (game.running) resetWaveAfterHit();
  }
}

/* ============================================================
   BUNKERS — destructible cover, built the same way as the invader
   sprites (a pixel grid), so damage is just removing cells.
   ============================================================ */
function initBunkers() {
  bunkers = [];
  var margin = 46;
  var spacing = (CANVAS_W - margin * 2) / (BUNKER_COUNT - 1);
  for (var i = 0; i < BUNKER_COUNT; i++) {
    var cells = BUNKER_TEMPLATE.map(function (row) {
      return row.split('').map(function (ch) { return ch === 'X'; });
    });
    bunkers.push({
      x: Math.round(margin + i * spacing - (BUNKER_COLS * BUNKER_PIXEL) / 2),
      y: Math.round(player.y - 96),
      cells: cells
    });
  }
}

// Removes a small cluster of cells around (x,y) if it lands on a live
// bunker cell. Returns true if it hit something (used to stop/consume
// non-piercing bullets from either side).
function damageBunkersAt(x, y) {
  var hit = false;
  bunkers.forEach(function (bk) {
    var col = Math.floor((x - bk.x) / BUNKER_PIXEL);
    var row = Math.floor((y - bk.y) / BUNKER_PIXEL);
    if (row >= 0 && row < BUNKER_ROWS && col >= 0 && col < BUNKER_COLS && bk.cells[row][col]) {
      hit = true;
      for (var dr = -1; dr <= 1; dr++) {
        for (var dc = -1; dc <= 1; dc++) {
          var rr = row + dr, cc = col + dc;
          if (rr >= 0 && rr < BUNKER_ROWS && cc >= 0 && cc < BUNKER_COLS) bk.cells[rr][cc] = false;
        }
      }
    }
  });
  return hit;
}

/* ============================================================
   PARTICLES — a short burst on invader kill, purely cosmetic.
   ============================================================ */
function spawnParticles(x, y, color) {
  for (var i = 0; i < 8; i++) {
    var ang = Math.random() * Math.PI * 2;
    var speed = 0.6 + Math.random() * 1.8;
    particles.push({
      x: x, y: y,
      vx: Math.cos(ang) * speed, vy: Math.sin(ang) * speed,
      life: 18 + Math.random() * 10, maxLife: 28,
      color: color
    });
  }
}

function updateParticles() {
  particles.forEach(function (p) {
    p.x += p.vx; p.y += p.vy;
    p.vx *= 0.94; p.vy *= 0.94;
    p.life--;
  });
  particles = particles.filter(function (p) { return p.life > 0; });
}

/* ============================================================
   SCREEN SHAKE / BANNERS / INCOMING CUES — small shared systems
   used by wave-clear, life-lost, and the capsule/shield warnings.
   ============================================================ */
function triggerShake() {
  shakeFramesLeft = SHAKE_FRAMES;
}

function startBanner(title, subtitle, frames, tone, onComplete) {
  game.banner = { title: title, subtitle: subtitle || '', framesLeft: frames, totalFrames: frames, tone: tone || 'wave', onComplete: onComplete || null, blocking: true };
}

function showIncomingCue(text, color) {
  incomingCue = { text: text, color: color || '#F5D76E', framesLeft: INCOMING_LEAD_FRAMES };
}

function updateShakeBannerAndCues() {
  if (shakeFramesLeft > 0) shakeFramesLeft--;

  if (game.banner) {
    game.banner.framesLeft--;
    if (game.banner.framesLeft <= 0) {
      var cb = game.banner.onComplete;
      game.banner = null;
      if (cb) cb();
    }
  }

  if (incomingCue) {
    incomingCue.framesLeft--;
    if (incomingCue.framesLeft <= 0) incomingCue = null;
  }

  updateParticles();
}

/* ============================================================
   PLAYER + BULLETS
   ============================================================ */
function updatePlayer() {
  if (keys.left) player.x -= PLAYER_SPEED;
  if (keys.right) player.x += PLAYER_SPEED;
  player.x = Math.max(20, Math.min(CANVAS_W - 20, player.x));
  if (fireCooldown > 0) fireCooldown--;
}

function tryFire() {
  if (!game.running || game.manualPaused || game.quizPaused || (game.banner && game.banner.blocking)) return;
  if (fireCooldown > 0) return;
  fireCooldown = FIRE_COOLDOWN_FRAMES;
  bullets.push({ x: player.x, y: player.y - 16, tier: game.bulletTier, pierce: game.bulletTier === 2, hitSet: {} });
  SFX.fire();
}

function updateBullets() {
  bullets.forEach(function (b) { b.y -= PLAYER_BULLET_SPEED; });
  bullets = bullets.filter(function (b) { return b.y > -20; });

  bullets.forEach(function (b) {
    var bunkerHit = damageBunkersAt(b.x, b.y);
    if (bunkerHit && !b.pierce) { b.dead = true; return; }

    invaders.forEach(function (inv, idx) {
      if (!inv.alive) return;
      if (b.pierce && b.hitSet[idx]) return;
      var p = invaderPos(inv);
      if (b.x > p.x - 2 && b.x < p.x + SPRITE_W + 2 && b.y > p.y - 2 && b.y < p.y + SPRITE_H + 2) {
        inv.alive = false;
        game.liveScore += (b.tier + 1) * 10;
        updateHud();
        spawnParticles(p.x + SPRITE_W / 2, p.y + SPRITE_H / 2, inv.row % 2 === 0 ? '#8C7CF0' : '#2FE0C4');
        SFX.pop();
        if (b.pierce) { b.hitSet[idx] = true; } else { b.dead = true; }
      }
    });
    if (capsule && capsule.active) {
      var dx = b.x - capsule.x, dy = b.y - capsule.y;
      if (dx * dx + dy * dy < 15 * 15) {
        capsule.active = false;
        capsuleState = 'idle'; capsuleTimer = 0;
        triggerQuestion();
      }
    }
    if (shieldBonus && shieldBonus.active) {
      var dx2 = b.x - shieldBonus.x, dy2 = b.y - shieldBonus.y;
      if (dx2 * dx2 + dy2 * dy2 < 14 * 14) {
        shieldBonus.active = false;
        shieldState = 'idle'; shieldTimer = 0;
        game.hasShield = true;
        SFX.shield();
        if (!b.pierce) b.dead = true;
      }
    }
  });
  bullets = bullets.filter(function (b) { return !b.dead; });
}

/* ============================================================
   ALIEN FIRE — only the frontmost alive invader per column can
   fire, at a low/occasional rate; shots are blockable by bunkers.
   ============================================================ */
function frontmostAliveByColumn() {
  var byCol = {};
  invaders.forEach(function (inv) {
    if (!inv.alive) return;
    if (!byCol[inv.col] || inv.row > byCol[inv.col].row) byCol[inv.col] = inv;
  });
  return Object.keys(byCol).map(function (k) { return byCol[k]; });
}

function updateAlienFire() {
  alienFireTimer++;
  if (alienFireTimer >= ALIEN_FIRE_CHECK_FRAMES) {
    alienFireTimer = 0;
    if (Math.random() < ALIEN_FIRE_CHANCE) {
      var candidates = frontmostAliveByColumn();
      if (candidates.length) {
        var shooter = candidates[Math.floor(Math.random() * candidates.length)];
        var p = invaderPos(shooter);
        invaderBullets.push({ x: p.x + SPRITE_W / 2, y: p.y + SPRITE_H });
      }
    }
  }

  invaderBullets.forEach(function (b) { b.y += ALIEN_BULLET_SPEED; });
  invaderBullets = invaderBullets.filter(function (b) { return b.y < CANVAS_H + 20; });

  invaderBullets.forEach(function (b) {
    if (damageBunkersAt(b.x, b.y)) { b.dead = true; return; }
    var dx = b.x - player.x, dy = b.y - player.y;
    if (dx * dx + dy * dy < 14 * 14) {
      b.dead = true;
      if (game.hasShield) {
        game.hasShield = false;
      } else {
        loseLife();
      }
    }
  });
  invaderBullets = invaderBullets.filter(function (b) { return !b.dead; });
}

/* ============================================================
   REAGENT CAPSULE
   ============================================================ */
function updateCapsule() {
  if (game.askedCount >= PUZZLES.length) return;

  if (capsuleState === 'active') {
    capsule.y += 0.7;
    if (capsule.y > player.y - 30) {
      capsule.active = false;
      capsuleState = 'idle'; capsuleTimer = 0;
      triggerQuestion();
    }
    return;
  }

  if (capsuleState === 'warning') {
    capsuleTimer--;
    if (capsuleTimer <= 0) {
      capsuleState = 'active';
      capsule = { x: 44 + Math.random() * (CANVAS_W - 88), y: 0, active: true };
    }
    return;
  }

  capsuleTimer++;
  if (capsuleTimer >= CAPSULE_INTERVAL_FRAMES - INCOMING_LEAD_FRAMES) {
    capsuleState = 'warning';
    capsuleTimer = INCOMING_LEAD_FRAMES;
    showIncomingCue('Reagent incoming', '#F5D76E');
  }
}

/* ============================================================
   SHIELD BONUS — a separate periodic pickup, unrelated to the quiz.
   Catch it (touch it with the ship, or shoot it) for one hit of
   protection against alien fire. Purely optional: if it isn't
   caught, it just drifts off the bottom with no penalty.
   ============================================================ */
function updateShieldBonus() {
  if (shieldState === 'active') {
    shieldBonus.y += SHIELD_FALL_SPEED;
    var dx = shieldBonus.x - player.x, dy = shieldBonus.y - player.y;
    if (dx * dx + dy * dy < 20 * 20) {
      shieldBonus.active = false;
      shieldState = 'idle'; shieldTimer = 0;
      game.hasShield = true;
      SFX.shield();
    } else if (shieldBonus.y > CANVAS_H + 16) {
      shieldBonus.active = false;
      shieldState = 'idle'; shieldTimer = 0;
    }
    return;
  }

  if (shieldState === 'warning') {
    shieldTimer--;
    if (shieldTimer <= 0) {
      shieldState = 'active';
      shieldBonus = { x: 40 + Math.random() * (CANVAS_W - 80), y: 0, active: true };
    }
    return;
  }

  shieldTimer++;
  if (shieldTimer >= SHIELD_INTERVAL_FRAMES - INCOMING_LEAD_FRAMES) {
    shieldState = 'warning';
    shieldTimer = INCOMING_LEAD_FRAMES;
    showIncomingCue('Shield incoming', SHIELD_COLOR);
  }
}

/* ============================================================
   QUIZ OVERLAY (one question per capsule, sequential, no repeats)
   ============================================================ */
function triggerQuestion() {
  if (activeQuestionIdx !== null || game.askedCount >= PUZZLES.length) return;
  game.quizPaused = true;
  var pauseBtn = $('btn-pause');
  if (pauseBtn) pauseBtn.classList.add('hidden');
  var pIdx = questionOrder[game.askedCount];
  activeQuestionIdx = pIdx;
  var puzzle = PUZZLES[pIdx];

  $('q-count').textContent = 'Question ' + (game.askedCount + 1) + ' of ' + PUZZLES.length;
  $('q-title').textContent = puzzle.title;
  $('q-reaction').innerHTML = renderReaction(puzzle.start, puzzle.target);

  var shuffled = shuffle(puzzle.options);
  var optsEl = $('q-options');
  optsEl.innerHTML = '';
  $('q-feedback').classList.add('hidden');

  shuffled.forEach(function (opt) {
    var btn = document.createElement('button');
    btn.className = 'option-card reagent-option';
    btn.disabled = true;
    btn.innerHTML = '<span class="reagent-text">' + opt.reagents + '</span>';
    btn.addEventListener('click', function () { answerQuestion(opt, btn, optsEl); });
    optsEl.appendChild(btn);
  });

  $('quiz-overlay').classList.remove('hidden');
  document.querySelector('.game-wrap').classList.add('hidden');
  document.querySelector('.controls').classList.add('hidden');

  // Guard against the tap/click that just shot the capsule bleeding into an
  // option button that renders in roughly the same screen spot a moment
  // later — options stay disabled for a beat after the overlay appears.
  setTimeout(function () {
    optsEl.querySelectorAll('.option-card').forEach(function (b) { b.disabled = false; });
  }, 450);
}

function answerQuestion(opt, btnEl, optsEl) {
  optsEl.querySelectorAll('.option-card').forEach(function (b) { b.disabled = true; });
  btnEl.classList.add(opt.correct ? 'correct' : 'incorrect');
  $('q-feedback-text').textContent = opt.explain;
  $('q-feedback').classList.remove('hidden');

  var puzzle = PUZZLES[activeQuestionIdx];
  var correctOpt = puzzle.options.filter(function (o) { return o.correct; })[0];
  game.picks[activeQuestionIdx] = { picked: opt.reagents, correct: correctOpt.reagents };

  if (opt.correct) {
    game.solved[activeQuestionIdx] = true;
    game.bulletTier = Math.min(2, game.bulletTier + 1);
    game.liveScore += 25;
    SFX.correct();
  } else {
    game.bulletTier = Math.max(0, game.bulletTier - 1);
    SFX.wrong();
  }
  updateHud();
}

function closeQuestionOverlay() {
  $('quiz-overlay').classList.add('hidden');
  document.querySelector('.game-wrap').classList.remove('hidden');
  document.querySelector('.controls').classList.remove('hidden');
  var pauseBtn = $('btn-pause');
  if (pauseBtn) pauseBtn.classList.remove('hidden');

  var wasCorrect = game.solved[activeQuestionIdx];
  game.askedCount++;
  activeQuestionIdx = null;
  game.quizPaused = false;
  updateHud();

  if (!wasCorrect) loseLife();
  if (game.running && game.askedCount >= PUZZLES.length) endGame(false);
}

/* ============================================================
   LIVES / END
   ============================================================ */
function loseLife() {
  game.lives--;
  updateHud();
  triggerShake();
  if (game.lives <= 0) { endGame(true); return; }
  SFX.lifeLost();
  startBanner('Life lost!', Math.max(game.lives, 0) + ' remaining', LIFE_BANNER_FRAMES, 'life');
}

function endGame(endedViaLivesLoss) {
  game.running = false;
  $('screen-game').classList.add('hidden');
  $('screen-end').classList.remove('hidden');

  var correctCount = game.solved.filter(Boolean).length;
  var pct = Math.round((correctCount / PUZZLES.length) * 100);
  var passed = !endedViaLivesLoss && pct >= 70;

  $('end-heading').textContent = passed ? 'Lab secured!' : (endedViaLivesLoss ? 'Out of lives.' : 'Worth another pass.');
  $('end-score').textContent = 'Questions correct: ' + correctCount + ' / ' + PUZZLES.length + ' (' + pct + '%). ' +
    'Lives remaining: ' + Math.max(game.lives, 0) + '. In-game score: ' + game.liveScore + ' (not sent to Moodle).';

  var reviewEl = $('end-review');
  reviewEl.innerHTML = '';
  for (var qi = 0; qi < game.askedCount; qi++) {
    var pIdx2 = questionOrder[qi];
    var wasSolved = game.solved[pIdx2];
    var pick = game.picks[pIdx2];
    var row = document.createElement('div');
    row.className = 'review-row ' + (wasSolved ? 'correct' : 'incorrect');
    var mainHtml = '<span class="review-mark">' + (wasSolved ? '\u2713' : '\u2717') + '</span><span>' + PUZZLES[pIdx2].title + '</span>';
    var detailHtml = '';
    if (!wasSolved && pick) {
      detailHtml = '<div class="review-detail">You chose: <strong>' + pick.picked + '</strong> \u2014 Correct: <strong>' + pick.correct + '</strong></div>';
    }
    row.innerHTML = '<div class="review-row-main">' + mainHtml + '</div>' + detailHtml;
    reviewEl.appendChild(row);
  }

  if (window.SCORM) {
    SCORM.setScore(pct, 100);
    SCORM.setComplete(passed);
  }
}

/* ============================================================
   HUD
   ============================================================ */
function updateHud() {
  $('hud-lives').textContent = '\u2665 '.repeat(Math.max(game.lives, 0)).trim() || 'No lives left';
  $('hud-tier').textContent = TIER_LABELS[game.bulletTier];
  $('hud-questions').textContent = 'Q ' + game.askedCount + ' / ' + PUZZLES.length;
  $('hud-score').textContent = 'Score ' + game.liveScore;
}

/* ============================================================
   RENDERING
   ============================================================ */
function initBackgroundDecor() {
  STAR_FIELD = [];
  for (var i = 0; i < 46; i++) {
    STAR_FIELD.push({ x: Math.random() * CANVAS_W, y: Math.random() * CANVAS_H, r: Math.random() * 1.3 + 0.4, tw: Math.random() * Math.PI * 2 });
  }
  BG_HEXES = [];
  for (var j = 0; j < 6; j++) {
    var r = 14 + Math.random() * 16;
    BG_HEXES.push({
      x: r + Math.random() * (CANVAS_W - 2 * r),
      y: r + Math.random() * (CANVAS_H - 2 * r),
      r: r
    });
  }
}

function drawBackground() {
  var grad = ctx.createLinearGradient(0, 0, 0, CANVAS_H);
  grad.addColorStop(0, '#1B1E45');
  grad.addColorStop(1, '#0A0C22');
  ctx.fillStyle = grad;
  ctx.fillRect(0, 0, CANVAS_W, CANVAS_H);

  ctx.save();
  ctx.globalAlpha = 0.07;
  ctx.strokeStyle = '#9AA0F5';
  ctx.lineWidth = 1.4;
  BG_HEXES.forEach(function (h) {
    ctx.beginPath();
    for (var k = 0; k < 6; k++) {
      var ang = (Math.PI / 180) * (60 * k - 90);
      var px = h.x + h.r * Math.cos(ang), py = h.y + h.r * Math.sin(ang);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.closePath();
    ctx.stroke();
  });
  ctx.restore();

  ctx.save();
  STAR_FIELD.forEach(function (s) {
    var tw = 0.5 + 0.5 * Math.sin(frameCount * 0.03 + s.tw);
    ctx.globalAlpha = 0.25 + 0.55 * tw;
    ctx.fillStyle = '#FFFFFF';
    ctx.beginPath();
    ctx.arc(s.x, s.y, s.r, 0, Math.PI * 2);
    ctx.fill();
  });
  ctx.restore();
}

function drawInvader(inv) {
  if (!inv.alive) return;
  var p = invaderPos(inv);
  var sprite = SPRITES[inv.variant];
  ctx.fillStyle = inv.row % 2 === 0 ? '#8C7CF0' : '#2FE0C4';
  for (var r = 0; r < sprite.length; r++) {
    for (var c = 0; c < sprite[r].length; c++) {
      if (sprite[r][c] === 'X') {
        ctx.fillRect(p.x + c * INV_PIXEL, p.y + r * INV_PIXEL, INV_PIXEL, INV_PIXEL);
      }
    }
  }
}

function drawBunkers() {
  bunkers.forEach(function (bk) {
    for (var r = 0; r < BUNKER_ROWS; r++) {
      for (var c = 0; c < BUNKER_COLS; c++) {
        if (bk.cells[r][c]) {
          ctx.fillStyle = '#B9C2F5';
          ctx.fillRect(bk.x + c * BUNKER_PIXEL, bk.y + r * BUNKER_PIXEL, BUNKER_PIXEL, BUNKER_PIXEL);
        }
      }
    }
  });
}

// A flickering thruster glow beneath the ship — brighter and larger while
// moving, a gentle idle pulse while still.
function drawEngineGlow(x, y) {
  var moving = keys.left || keys.right;
  var flicker = moving
    ? (0.6 + 0.3 * Math.sin(frameCount * 0.8) + 0.2 * Math.random())
    : (0.32 + 0.08 * Math.sin(frameCount * 0.15));
  var len = moving ? (11 + 4 * Math.random()) : 6;
  ctx.save();
  ctx.globalAlpha = Math.max(0.2, Math.min(1, flicker));
  var grad = ctx.createRadialGradient(x, y + 20, 1, x, y + 20 + len, len + 4);
  grad.addColorStop(0, '#FFE9A8');
  grad.addColorStop(1, 'rgba(255, 107, 74, 0)');
  ctx.fillStyle = grad;
  ctx.beginPath();
  ctx.ellipse(x, y + 20 + len / 2, 5, len, 0, 0, Math.PI * 2);
  ctx.fill();
  ctx.restore();
}

// A small flask riding a rocket frame — reads as both a ship (nose cone,
// delta wings, engine glow) and a piece of glassware (neck, flared body,
// rounded base). The "liquid" fill tracks the current bullet tier, so the
// ship itself shows the player's upgrade state at a glance.
function drawPlayer() {
  var x = player.x, y = player.y;

  drawEngineGlow(x, y);

  ctx.fillStyle = 'rgba(245, 215, 110, 0.5)';
  ctx.beginPath();
  ctx.ellipse(x, y + 15, 6.5, 9, 0, 0, Math.PI * 2);
  ctx.fill();

  ctx.beginPath();
  ctx.moveTo(x - 4, y - 12);
  ctx.lineTo(x + 4, y - 12);
  ctx.lineTo(x + 4, y - 5);
  ctx.lineTo(x + 13, y + 9);
  ctx.quadraticCurveTo(x, y + 15, x - 13, y + 9);
  ctx.lineTo(x - 4, y - 5);
  ctx.closePath();
  ctx.fillStyle = '#E9E7F7';
  ctx.fill();

  ctx.save();
  ctx.beginPath();
  ctx.moveTo(x - 4, y - 2);
  ctx.lineTo(x + 4, y - 2);
  ctx.lineTo(x + 13, y + 9);
  ctx.quadraticCurveTo(x, y + 15, x - 13, y + 9);
  ctx.closePath();
  ctx.clip();
  ctx.fillStyle = TIER_COLORS[game.bulletTier];
  ctx.fillRect(x - 16, y - 4, 32, 22);
  ctx.restore();

  ctx.strokeStyle = '#B7B3D9';
  ctx.lineWidth = 1.2;
  ctx.beginPath();
  ctx.moveTo(x - 4, y - 12);
  ctx.lineTo(x + 4, y - 12);
  ctx.lineTo(x + 4, y - 5);
  ctx.lineTo(x + 13, y + 9);
  ctx.quadraticCurveTo(x, y + 15, x - 13, y + 9);
  ctx.lineTo(x - 4, y - 5);
  ctx.closePath();
  ctx.stroke();

  ctx.fillStyle = '#FF6B4A';
  ctx.beginPath();
  ctx.moveTo(x - 5, y - 12);
  ctx.lineTo(x + 5, y - 12);
  ctx.lineTo(x, y - 20);
  ctx.closePath();
  ctx.fill();

  ctx.fillStyle = '#4A38C4';
  ctx.beginPath();
  ctx.moveTo(x - 10, y + 4);
  ctx.lineTo(x - 19, y + 12);
  ctx.lineTo(x - 8, y + 10);
  ctx.closePath();
  ctx.fill();
  ctx.beginPath();
  ctx.moveTo(x + 10, y + 4);
  ctx.lineTo(x + 19, y + 12);
  ctx.lineTo(x + 8, y + 10);
  ctx.closePath();
  ctx.fill();
}

function drawBullet(b) {
  ctx.strokeStyle = TIER_COLORS[b.tier];
  ctx.lineWidth = 2;
  var lines = b.tier + 1;
  var spacing = 4;
  var startOffset = -((lines - 1) * spacing) / 2;
  for (var i = 0; i < lines; i++) {
    var lx = b.x + startOffset + i * spacing;
    ctx.beginPath();
    ctx.moveTo(lx, b.y - 8);
    ctx.lineTo(lx, b.y + 8);
    ctx.stroke();
  }
  // Each bullet is the simplest real molecule with that bond order —
  // ethane / ethene / ethyne — with its substituent labelled at both ends.
  ctx.fillStyle = TIER_COLORS[b.tier];
  ctx.font = 'bold 8px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText(BULLET_LABELS[b.tier], b.x, b.y - 9);
  ctx.textBaseline = 'top';
  ctx.fillText(BULLET_LABELS[b.tier], b.x, b.y + 9);
}

function drawAlienBullet(b) {
  ctx.strokeStyle = '#FF6B4A';
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(b.x, b.y - 6);
  ctx.lineTo(b.x, b.y + 6);
  ctx.stroke();
}

function drawCapsule() {
  if (!capsule || !capsule.active) return;
  ctx.fillStyle = '#F5D76E';
  ctx.beginPath();
  ctx.arc(capsule.x, capsule.y, 11, 0, Math.PI * 2);
  ctx.fill();
  ctx.fillStyle = '#10132A';
  ctx.font = 'bold 13px sans-serif';
  ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillText('?', capsule.x, capsule.y + 1);
}

function hexPath(x, y, r) {
  ctx.beginPath();
  for (var k = 0; k < 6; k++) {
    var ang = (Math.PI / 180) * (60 * k - 90);
    var px = x + r * Math.cos(ang), py = y + r * Math.sin(ang);
    if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
  }
  ctx.closePath();
}

function shieldGlyphPath(x, y, s) {
  ctx.beginPath();
  ctx.moveTo(x, y - s);
  ctx.lineTo(x + s * 0.8, y - s * 0.55);
  ctx.lineTo(x + s * 0.8, y + s * 0.25);
  ctx.quadraticCurveTo(x + s * 0.8, y + s, x, y + s * 1.25);
  ctx.quadraticCurveTo(x - s * 0.8, y + s, x - s * 0.8, y + s * 0.25);
  ctx.lineTo(x - s * 0.8, y - s * 0.55);
  ctx.closePath();
}

// A shield icon inside a small hexagon — a separate, purely optional
// pickup (not part of the quiz), caught by touching it with the ship or
// shooting it.
function drawShieldBonus() {
  if (!shieldBonus || !shieldBonus.active) return;
  ctx.save();
  hexPath(shieldBonus.x, shieldBonus.y, 12);
  ctx.fillStyle = 'rgba(79, 216, 232, 0.18)';
  ctx.fill();
  ctx.strokeStyle = SHIELD_COLOR;
  ctx.lineWidth = 2;
  ctx.stroke();
  shieldGlyphPath(shieldBonus.x, shieldBonus.y, 5.5);
  ctx.fillStyle = SHIELD_COLOR;
  ctx.fill();
  ctx.restore();
}

// The active-shield indicator around the ship — same hex language as the
// pickup, so it's clear at a glance what protection is currently up.
function drawShieldRing() {
  if (!game.hasShield) return;
  ctx.save();
  var pulse = 1 + 0.06 * Math.sin(frameCount * 0.12);
  hexPath(player.x, player.y - 2, 24 * pulse);
  ctx.strokeStyle = SHIELD_COLOR;
  ctx.lineWidth = 2;
  ctx.globalAlpha = 0.85;
  ctx.stroke();
  ctx.fillStyle = 'rgba(79, 216, 232, 0.08)';
  ctx.fill();
  ctx.restore();
}

function drawParticles() {
  particles.forEach(function (p) {
    var alpha = Math.max(0, p.life / p.maxLife);
    ctx.globalAlpha = alpha;
    ctx.fillStyle = p.color;
    ctx.fillRect(p.x - 1.5, p.y - 1.5, 3, 3);
  });
  ctx.globalAlpha = 1;
}

function roundRectPath(x, y, w, h, r) {
  ctx.beginPath();
  ctx.moveTo(x + r, y);
  ctx.arcTo(x + w, y, x + w, y + h, r);
  ctx.arcTo(x + w, y + h, x, y + h, r);
  ctx.arcTo(x, y + h, x, y, r);
  ctx.arcTo(x, y, x + w, y, r);
  ctx.closePath();
}

// A brief heads-up pill near the top of the canvas before a capsule or
// shield drops — non-blocking, gameplay keeps moving underneath it.
function drawIncomingCue() {
  if (!incomingCue) return;
  var alpha = Math.min(1, incomingCue.framesLeft / 15);
  ctx.save();
  ctx.globalAlpha = alpha;
  ctx.font = 'bold 12px sans-serif';
  var text = incomingCue.text;
  var w = ctx.measureText(text).width + 24;
  var cx = CANVAS_W / 2, cy = 26;
  ctx.fillStyle = 'rgba(10, 12, 34, 0.75)';
  roundRectPath(cx - w / 2, cy - 12, w, 24, 12);
  ctx.fill();
  ctx.strokeStyle = incomingCue.color;
  ctx.lineWidth = 1.5;
  ctx.stroke();
  ctx.fillStyle = incomingCue.color;
  ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillText(text, cx, cy + 1);
  ctx.restore();
}

// The wave-cleared / life-lost banner — a brief centred callout with a
// soft fade in and out, drawn outside the shake transform so its text
// stays legible even while the world around it is shaking.
function drawBanner() {
  if (!game.banner) return;
  var b = game.banner;
  var t = 1 - b.framesLeft / b.totalFrames;
  var alpha = Math.min(1, Math.min(t, 1 - t) * 6);
  ctx.save();
  ctx.globalAlpha = alpha * 0.55;
  ctx.fillStyle = '#05061A';
  ctx.fillRect(0, CANVAS_H / 2 - 46, CANVAS_W, 92);
  ctx.restore();

  ctx.save();
  ctx.globalAlpha = alpha;
  ctx.textAlign = 'center';
  ctx.fillStyle = b.tone === 'life' ? '#FF8A65' : '#8C7CF0';
  ctx.font = 'bold 20px sans-serif';
  ctx.fillText(b.title, CANVAS_W / 2, CANVAS_H / 2 - 6);
  if (b.subtitle) {
    ctx.fillStyle = '#CFD3F5';
    ctx.font = '13px sans-serif';
    ctx.fillText(b.subtitle, CANVAS_W / 2, CANVAS_H / 2 + 18);
  }
  ctx.restore();
}

function draw() {
  frameCount++;
  ctx.save();
  if (shakeFramesLeft > 0) {
    var mag = SHAKE_MAGNITUDE * (shakeFramesLeft / SHAKE_FRAMES);
    ctx.translate((Math.random() * 2 - 1) * mag, (Math.random() * 2 - 1) * mag);
  }
  drawBackground();
  drawBunkers();
  invaders.forEach(drawInvader);
  bullets.forEach(drawBullet);
  invaderBullets.forEach(drawAlienBullet);
  drawParticles();
  drawCapsule();
  drawShieldBonus();
  drawShieldRing();
  drawPlayer();
  ctx.restore();
  drawIncomingCue();
  drawBanner();
}

/* ============================================================
   MAIN LOOP
   ============================================================ */
function update() {
  // Manual pause freezes everything, including banner/shake timers, so
  // nothing progresses while the pause card is up.
  if (game.manualPaused) return;

  updateShakeBannerAndCues();
  if (game.quizPaused || (game.banner && game.banner.blocking)) return;

  updatePlayer();
  updateBullets();
  updateInvaders();
  if (!game.running) return;
  updateAlienFire();
  if (!game.running) return;
  updateCapsule();
  if (!game.running) return;
  updateShieldBonus();
}

function loop() {
  if (!game.running) return;
  update();
  draw();
  requestAnimationFrame(loop);
}

/* ============================================================
   START / RESTART
   ============================================================ */
function startGame() {
  SFX.unlock(); // Start is a user gesture — safe point to spin up the AudioContext

  var preset = DIFFICULTY_PRESETS[difficulty] || DIFFICULTY_PRESETS.normal;
  ALIEN_FIRE_CHANCE = preset.fireChance;
  INV_WAVE_STEP = preset.waveStep;

  game.lives = preset.lives;
  game.liveScore = 0;
  game.solved = [false, false, false, false, false, false, false, false];
  game.picks = [];
  game.askedCount = 0;
  game.bulletTier = 0;
  game.wave = 1;
  game.quizPaused = false;
  game.manualPaused = false;
  game.banner = null;
  game.running = true;

  questionOrder = shuffle(PUZZLES.map(function (_, i) { return i; }));
  activeQuestionIdx = null;
  bullets = [];
  invaderBullets = [];
  particles = [];
  incomingCue = null;
  shakeFramesLeft = 0;
  alienFireTimer = 0;
  frameCount = 0;
  capsule = null;
  capsuleState = 'idle'; capsuleTimer = 0;
  shieldBonus = null;
  shieldState = 'idle'; shieldTimer = 0;
  game.hasShield = false;
  fireCooldown = 0;
  player = { x: CANVAS_W / 2, y: CANVAS_H - 40 };
  initInvaders();
  initBunkers();
  initBackgroundDecor();

  $('screen-intro').classList.add('hidden');
  $('screen-end').classList.add('hidden');
  $('screen-game').classList.remove('hidden');
  $('quiz-overlay').classList.add('hidden');
  $('pause-overlay').classList.add('hidden');
  $('btn-pause').classList.remove('hidden');
  $('btn-pause').textContent = 'Pause';
  document.querySelector('.game-wrap').classList.remove('hidden');
  document.querySelector('.controls').classList.remove('hidden');

  updateHud();
  requestAnimationFrame(loop);
}

function togglePause() {
  if (!game.running || game.quizPaused) return;
  game.manualPaused = !game.manualPaused;
  $('pause-overlay').classList.toggle('hidden', !game.manualPaused);
  $('btn-pause').textContent = game.manualPaused ? 'Resume' : 'Pause';
}

/* ============================================================
   INPUT
   ============================================================ */
document.addEventListener('DOMContentLoaded', function () {
  canvas = $('game-canvas');
  ctx = canvas.getContext('2d');

  if (window.SCORM) SCORM.init();

  $('btn-start').addEventListener('click', startGame);
  $('btn-restart').addEventListener('click', startGame);
  $('btn-q-continue').addEventListener('click', closeQuestionOverlay);

  // Difficulty selector on the intro screen — Normal is the default.
  var diffButtons = document.querySelectorAll('.diff-btn');
  diffButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      difficulty = btn.getAttribute('data-diff');
      diffButtons.forEach(function (b) { b.classList.toggle('active', b === btn); });
    });
  });

  // Sound is off by default; the mute button just toggles SFX's own flag.
  var muteBtn = $('btn-mute');
  muteBtn.textContent = '♪ Off';
  muteBtn.setAttribute('aria-pressed', 'false');
  muteBtn.addEventListener('click', function () {
    var newMuted = !SFX.isMuted();
    SFX.setMuted(newMuted);
    muteBtn.textContent = newMuted ? '♪ Off' : '♪ On';
    muteBtn.setAttribute('aria-pressed', String(!newMuted));
  });

  $('btn-pause').addEventListener('click', togglePause);
  $('btn-resume').addEventListener('click', togglePause);

  document.addEventListener('keydown', function (e) {
    if (e.key === 'ArrowLeft' || e.key === 'a' || e.key === 'A') { keys.left = true; e.preventDefault(); }
    if (e.key === 'ArrowRight' || e.key === 'd' || e.key === 'D') { keys.right = true; e.preventDefault(); }
    if (e.key === ' ') { tryFire(); e.preventDefault(); }
    if (e.key === 'p' || e.key === 'P') { togglePause(); }
  });
  document.addEventListener('keyup', function (e) {
    if (e.key === 'ArrowLeft' || e.key === 'a' || e.key === 'A') keys.left = false;
    if (e.key === 'ArrowRight' || e.key === 'd' || e.key === 'D') keys.right = false;
  });

  var lb = $('ctl-left'), rb = $('ctl-right'), fb = $('ctl-fire');
  lb.addEventListener('touchstart', function (e) { keys.left = true; e.preventDefault(); }, { passive: false });
  lb.addEventListener('touchend', function (e) { keys.left = false; e.preventDefault(); }, { passive: false });
  lb.addEventListener('mousedown', function () { keys.left = true; });
  lb.addEventListener('mouseup', function () { keys.left = false; });
  rb.addEventListener('touchstart', function (e) { keys.right = true; e.preventDefault(); }, { passive: false });
  rb.addEventListener('touchend', function (e) { keys.right = false; e.preventDefault(); }, { passive: false });
  rb.addEventListener('mousedown', function () { keys.right = true; });
  rb.addEventListener('mouseup', function () { keys.right = false; });
  fb.addEventListener('touchstart', function (e) { tryFire(); e.preventDefault(); }, { passive: false });
  fb.addEventListener('mousedown', function () { tryFire(); });

  window.addEventListener('beforeunload', function () {
    if (window.SCORM) SCORM.finish();
  });
});
