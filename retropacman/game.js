/* ============================================================
   MOLECULE RENDERER  (v2 — true skeletal style)
   A molecule is { atoms: [...], bonds: [...] }.
   atoms[i] is one of:
     RING()            benzene ring icon
     C()               implicit carbon vertex (no label, plain kink)
     C(BR(...))         implicit carbon vertex carrying a branch group
     LBL('OH')         an explicit heteroatom / group, inline or terminal
   bonds[i] is 'single' or 'double' for the bond between atoms[i] and atoms[i+1]
   (defaults to 'single' if omitted).
   Every molecule is drawn at the SAME bond length / font size in real
   pixels, so target and option diagrams are always visually consistent.
   ============================================================ */

var BOND = 30;      // px, default backbone bond length (stretched if a wide label needs more room)
var RING_R = 15;     // px, ring icon radius
var BRANCH_LEN = 24; // px, default vertical branch bond length
var FONT_SIZE = 13;  // px, atom label font size (matches CSS .atom-label)
var DBL_OFFSET = 2.6; // px, gap between parallel lines of a double bond
var MIN_SEG = 7;      // px, minimum visible line length once end-padding is subtracted

function RING() { return { t: 'ring' }; }
function C(branch) { return branch ? { t: 'c', branch: branch } : { t: 'c' }; }
function LBL(text) { return { t: 'label', text: text }; }
function BR(text, dir, bond) { return { text: text, dir: dir || 'up', bond: bond || 'single' }; }
// C3: a standalone carbon fully drawn out with several substituents placed at
// explicit angles (degrees, 0=right, 90=down, -90=up), each 120 degrees apart —
// used for small molecules like formaldehyde where every substituent (incl. H)
// needs to be shown, rather than an implicit chain vertex.
function C3(branches) { return { t: 'c3', branches: branches }; }
function BR3(text, angleDeg, bond) { return { text: text, angle: angleDeg, bond: bond || 'single' }; }
function MOL(atoms, bonds) { return { atoms: atoms, bonds: bonds || [] }; }

function textHalfWidth(text) {
  return Math.max(6, text.length * 3.35 + 2);
}

function ringIconSVG(cx, cy, r) {
  var pts = [];
  for (var k = 0; k < 6; k++) {
    var ang = (Math.PI / 180) * (60 * k - 90);
    pts.push((cx + r * Math.cos(ang)).toFixed(1) + ',' + (cy + r * Math.sin(ang)).toFixed(1));
  }
  return '<polygon points="' + pts.join(' ') + '" class="ring-hex"/>' +
         '<circle cx="' + cx.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="' + (r * 0.55).toFixed(1) + '" class="ring-circle"/>';
}

function lineSVG(x1, y1, x2, y2, dashed) {
  return '<line x1="' + x1.toFixed(1) + '" y1="' + y1.toFixed(1) + '" x2="' + x2.toFixed(1) + '" y2="' + y2.toFixed(1) +
    '" stroke="' + (dashed ? 'var(--break)' : 'var(--ink)') + '" stroke-width="' + (dashed ? 2.4 : 2) +
    '" stroke-dasharray="' + (dashed ? '5,4' : '0') + '" stroke-linecap="round"/>';
}

function formatChemText(text) {
  var out = '';
  for (var i = 0; i < text.length; i++) {
    var ch = text[i];
    if (ch >= '0' && ch <= '9') {
      out += '<tspan baseline-shift="sub" font-size="72%">' + ch + '</tspan>';
    } else {
      out += ch;
    }
  }
  return out;
}

function labelSVG(x, y, text) {
  return '<text x="' + x.toFixed(1) + '" y="' + (y + FONT_SIZE * 0.36).toFixed(1) + '" text-anchor="middle" class="atom-label">' + formatChemText(text) + '</text>';
}

// End-padding (how far from the vertex centre a bond line should stop)
function endPad(atom) {
  if (atom.t === 'ring') return RING_R;
  if (atom.t === 'label') return textHalfWidth(atom.text);
  return 0; // plain implicit-carbon vertex — line runs straight through it
}

// Determines, for each backbone bond, whether it steps "up" (-1) or "down" (+1).
// A branched atom (carbonyl, OH, etc.) must sit at a true peak/valley so its two
// backbone neighbours both angle away from the branch at 120 degrees (branch up
// -> both neighbours go down; branch down -> both neighbours go up).
function computeSegSigns(atoms, bonds) {
  var n = atoms.length;
  var branchIdx = -1;
  for (var i = 0; i < n; i++) {
    if (atoms[i].t === 'c' && atoms[i].branch) { branchIdx = i; break; }
  }
  var s0 = 1;
  if (branchIdx > 0) {
    var required = (atoms[branchIdx].branch.dir === 'down') ? 1 : -1; // sign of segSign[branchIdx-1]
    var parityPow = ((branchIdx - 1) % 2 === 0) ? 1 : -1;
    s0 = required * parityPow;
  }
  var segSign = [];
  for (i = 0; i < n - 1; i++) segSign.push(s0 * ((i % 2 === 0) ? 1 : -1));

  // A triple bond is linear (sp carbon, 180 degrees) — the bond itself and
  // both its neighbouring bonds must be collinear. Rather than forcing them
  // flat (which would land a ring attachment in the gap between two vertices,
  // since the hexagon's vertices sit at +-30/+-90/+-150 degrees, not 0), they
  // all take the SAME sign as the bond immediately before the triple bond, so
  // the whole straight run continues at the usual +-30 degree bond angle —
  // still perfectly collinear, but landing exactly on a ring vertex.
  bonds = bonds || [];
  for (i = 0; i < n - 1; i++) {
    if (bonds[i] === 'triple') {
      var s = (i - 1 >= 0) ? segSign[i - 1] : segSign[i];
      segSign[i] = s;
      if (i + 1 <= n - 2) segSign[i + 1] = s;
    }
  }
  return segSign;
}

function drawBond(x1, y1, x2, y2, pad1, pad2, order, dashed) {
  order = order === true ? 'double' : (order === false ? 'single' : (order || 'single'));
  var dx = x2 - x1, dy = y2 - y1;
  var len = Math.sqrt(dx * dx + dy * dy) || 1;
  var ux = dx / len, uy = dy / len;
  // Guarantee a minimum visible segment even if padding is large relative to length
  if (pad1 + pad2 + MIN_SEG > len) {
    var scale = (len - MIN_SEG) / (pad1 + pad2 || 1);
    if (scale < 0) scale = 0;
    pad1 *= scale; pad2 *= scale;
  }
  var ax = x1 + ux * pad1, ay = y1 + uy * pad1;
  var bx = x2 - ux * pad2, by = y2 - uy * pad2;
  if (order === 'single') return lineSVG(ax, ay, bx, by, dashed);
  if (order === 'double') {
    var px = -uy * DBL_OFFSET, py = ux * DBL_OFFSET;
    return lineSVG(ax + px, ay + py, bx + px, by + py, dashed) +
           lineSVG(ax - px, ay - py, bx - px, by - py, dashed);
  }
  // triple: centre line plus two parallel outer lines
  var px2 = -uy * DBL_OFFSET * 1.15, py2 = ux * DBL_OFFSET * 1.15;
  return lineSVG(ax, ay, bx, by, dashed) +
         lineSVG(ax + px2, ay + py2, bx + px2, by + py2, dashed) +
         lineSVG(ax - px2, ay - py2, bx - px2, by - py2, dashed);
}

function renderMol(mol, opts) {
  opts = opts || {};
  var atoms = mol.atoms, bonds = mol.bonds || [];
  var n = atoms.length;
  var highlightBond = opts.highlightBond;
  var segSign = computeSegSigns(atoms, bonds);

  // Each bond gets its own length: the default BOND, or longer if the two
  // atoms it connects need more clearance for their labels (avoids the
  // ring/label overlap and vanishing-line bugs with wide labels).
  var pos = [{ x: 0, y: 0 }];
  for (var i = 0; i < n - 1; i++) {
    var pad1 = endPad(atoms[i]), pad2 = endPad(atoms[i + 1]);
    var segLen = Math.max(BOND, pad1 + pad2 + MIN_SEG + 4);
    var sdx, sdy;
    if (segSign[i] === 0) {
      // part of a linear (triple-bond) run — straight, full-length segment
      sdx = segLen; sdy = 0;
    } else {
      sdx = segLen * Math.cos(Math.PI / 6);
      sdy = segSign[i] * segLen * Math.sin(Math.PI / 6);
    }
    pos.push({ x: pos[i].x + sdx, y: pos[i].y + sdy });
  }

  var body = '';
  var minX = 0, maxX = 0, minY = 0, maxY = 0;

  function track(x, y) {
    if (x < minX) minX = x; if (x > maxX) maxX = x;
    if (y < minY) minY = y; if (y > maxY) maxY = y;
  }
  for (i = 0; i < n; i++) {
    var pad = endPad(atoms[i]);
    track(pos[i].x - pad, pos[i].y - pad);
    track(pos[i].x + pad, pos[i].y + pad);
  }

  // backbone bonds
  for (i = 0; i < n - 1; i++) {
    var order = bonds[i] || 'single';
    var isHi = (i === highlightBond);
    body += drawBond(pos[i].x, pos[i].y, pos[i + 1].x, pos[i + 1].y, endPad(atoms[i]), endPad(atoms[i + 1]), order, isHi);
  }

  // branches (single substituent hanging off a backbone carbon, e.g. C=O, OH)
  for (i = 0; i < n; i++) {
    var a = atoms[i];
    if (a.t === 'c' && a.branch) {
      var vert = (a.branch.dir === 'down') ? 1 : -1;
      var labelPad = textHalfWidth(a.branch.text);
      var branchLen = Math.max(BRANCH_LEN, labelPad + MIN_SEG + 4);
      var bx = pos[i].x, by = pos[i].y + vert * branchLen;
      body += drawBond(pos[i].x, pos[i].y, bx, by, 0, labelPad, a.branch.bond || 'single', false);
      body += labelSVG(bx, by, a.branch.text);
      track(bx - labelPad, by - 10);
      track(bx + labelPad, by + 10);
    }
    // c3: a standalone, fully-drawn-out carbon (e.g. formaldehyde) with several
    // substituents at explicit angles, each 120 degrees apart from the next.
    if (a.t === 'c3') {
      a.branches.forEach(function (br) {
        var rad = br.angle * Math.PI / 180;
        var lp = textHalfWidth(br.text);
        var bl = Math.max(BRANCH_LEN, lp + MIN_SEG + 4);
        var ex = pos[i].x + bl * Math.cos(rad);
        var ey = pos[i].y + bl * Math.sin(rad);
        body += drawBond(pos[i].x, pos[i].y, ex, ey, 0, lp, br.bond || 'single', false);
        body += labelSVG(ex, ey, br.text);
        track(ex - lp, ey - 10);
        track(ex + lp, ey + 10);
      });
    }
  }

  // atom glyphs (rings + labels) drawn after bonds so they sit on top
  for (i = 0; i < n; i++) {
    if (atoms[i].t === 'ring') {
      body += ringIconSVG(pos[i].x, pos[i].y, RING_R);
    } else if (atoms[i].t === 'label') {
      body += labelSVG(pos[i].x, pos[i].y, atoms[i].text);
    }
  }

  var margin = 8;
  var w = (maxX - minX) + margin * 2;
  var h = (maxY - minY) + margin * 2;
  var offX = margin - minX, offY = margin - minY;

  return '<svg viewBox="0 0 ' + w.toFixed(1) + ' ' + h.toFixed(1) + '" width="' + w.toFixed(0) + '" height="' + h.toFixed(0) +
    '" xmlns="http://www.w3.org/2000/svg" class="mol"><g transform="translate(' + offX.toFixed(1) + ',' + offY.toFixed(1) + ')">' +
    body + '</g></svg>';
}

function renderPrecursorPair(mol1, mol2) {
  var html = '<div class="option-precursors">';
  html += renderMol(mol1);
  if (mol2) {
    html += '<span class="plus-sign">+</span>';
    html += renderMol(mol2);
  }
  html += '</div>';
  return html;
}

/* ============================================================
   PUZZLE DATA — all 8, split 4/4 across two linked maze rooms
   ============================================================ */

var PUZZLES = [
  {
    title: 'Ethyl benzoate',
    target: MOL([RING(), C(BR('O', 'up', 'double')), LBL('O'), C(), C()]),
    highlightBond: 1,
    options: [
      {
        correct: true,
        mol1: MOL([RING(), C(BR('O', 'up', 'double')), LBL('OH')]),
        mol2: MOL([C(), C(), LBL('OH')]),
        condition: 'H⁺, heat',
        label: 'Benzoic acid + Ethanol',
        explain: 'Correct — disconnecting the acyl–oxygen bond gives benzoic acid and ethanol. Forward direction: Fischer esterification (acid + alcohol, H⁺ catalyst).'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(), LBL('OH')]),
        mol2: MOL([C(), C(BR('O', 'up', 'double')), LBL('OH')]),
        label: 'Benzyl alcohol + Acetic acid',
        explain: 'The acyl and alkyl portions are swapped — this pair would give benzyl acetate, not ethyl benzoate.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('OH')]),
        mol2: MOL([C(), C(), C(BR('O', 'up', 'double')), LBL('OH')]),
        label: 'Phenol + Propanoic acid',
        explain: 'Wrong connectivity and chain length — this would give phenyl propanoate, a different ester entirely.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(BR('O', 'up', 'double')), LBL('Cl')]),
        mol2: MOL([C(), C(), C(), LBL('OH')]),
        label: 'Benzoyl chloride + Propan-1-ol',
        explain: 'The acyl chloride is a reasonable acylating agent, but the alcohol has one carbon too many — this gives propyl benzoate.'
      }
    ]
  },
  {
    title: 'Acetanilide',
    target: MOL([RING(), LBL('NH'), C(BR('O', 'up', 'double')), C()]),
    highlightBond: 1,
    options: [
      {
        correct: true,
        mol1: MOL([RING(), LBL('NH2')]),
        mol2: MOL([C(), C(BR('O', 'up', 'double')), LBL('Cl')]),
        condition: 'base',
        label: 'Aniline + Acetyl chloride',
        explain: 'Correct — disconnecting the N–C(=O) bond gives aniline and acetyl chloride. Forward direction: acylation of the amine.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(BR('O', 'up', 'double')), LBL('NH2')]),
        mol2: MOL([C(), LBL('NH2')]),
        label: 'Benzamide + Methylamine',
        explain: 'This breaks a different bond and swaps which fragment carries the carbonyl — it leads to N-methylbenzamide, an isomeric but different amide.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('NO2')]),
        mol2: MOL([C(), C(BR('O', 'up', 'double')), LBL('Cl')]),
        label: 'Nitrobenzene + Acetyl chloride',
        explain: 'Nitrobenzene has the wrong oxidation state at nitrogen — it would need reduction to aniline first before this disconnection applies.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('OH')]),
        mol2: MOL([C(), C(BR('O', 'up', 'double')), LBL('NH2')]),
        label: 'Phenol + Acetamide',
        explain: 'This pair of functional groups does not combine to form the marked N–C bond — it points toward an ester, not an amide.'
      }
    ]
  },
  {
    title: '1-Phenylethanol',
    target: MOL([RING(), C(BR('OH', 'up', 'single')), C()]),
    highlightBond: 1,
    options: [
      {
        correct: true,
        mol1: MOL([RING(), C(BR('O', 'up', 'double'))]),
        mol2: MOL([C(), LBL('MgBr')]),
        condition: 'then H₃O⁺',
        label: 'Benzaldehyde + Methylmagnesium bromide',
        explain: 'Correct — the C–CH3 bond next to the alcohol is disconnected as a Grignard addition: the methyl nucleophile adds to the aldehyde carbonyl.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(BR('O', 'up', 'double')), C()]),
        label: 'Acetophenone (alone)',
        explain: 'This is the oxidised ketone, not a disconnection — reducing it would work synthetically, but it does not break the marked C–C bond retrosynthetically.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('MgBr')]),
        mol2: MOL([C(), C(BR('O', 'up', 'double'))]),
        label: 'Phenylmagnesium bromide + Acetaldehyde',
        explain: 'This breaks the Ph–C bond instead of the marked C–CH3 bond — a valid Grignard route to the same product, but not the disconnection shown here.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(), C()]),
        label: 'Ethylbenzene (alone)',
        explain: 'There is no bond-forming logic here — this is simply the deoxygenated hydrocarbon, not a retrosynthetic precursor pair.'
      }
    ]
  },
  {
    title: 'Diphenylmethanol',
    target: MOL([RING(), C(BR('OH', 'up', 'single')), RING()]),
    highlightBond: 1,
    options: [
      {
        correct: true,
        mol1: MOL([RING(), C(BR('O', 'up', 'double'))]),
        mol2: MOL([RING(), LBL('MgBr')]),
        condition: 'then H₃O⁺',
        label: 'Benzaldehyde + Phenylmagnesium bromide',
        explain: 'Correct — disconnecting the marked C–Ph bond gives benzaldehyde and a phenyl Grignard, which adds to the carbonyl.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(BR('O', 'up', 'double')), RING()]),
        label: 'Benzophenone (alone)',
        explain: 'This is the oxidised ketone. Reduction (e.g. NaBH4) gets you to the product, but it is not a disconnection of the marked bond.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(), RING()]),
        label: 'Diphenylmethane (alone)',
        explain: 'This lacks the oxygen entirely and has no bond-forming step that reconnects to give the alcohol.'
      },
      {
        correct: false,
        mol1: MOL([RING(), RING()]),
        mol2: MOL([C3([BR3('O', -90, 'double'), BR3('H', 30, 'single'), BR3('H', 150, 'single')])]),
        label: 'Biphenyl + Formaldehyde',
        explain: 'Wrong connectivity — biphenyl already has the two rings joined directly, which is not the bond pattern in the target.'
      }
    ]
  },
  {
    title: '(E)-Stilbene',
    target: MOL([RING(), C(), C(), RING()], ['single', 'double', 'single']),
    highlightBond: 1,
    options: [
      {
        correct: true,
        mol1: MOL([RING(), C(BR('O', 'up', 'double'))]),
        mol2: MOL([RING(), C(BR('PPh3', 'up', 'double'))]),
        condition: 'Wittig',
        label: 'Benzaldehyde + a benzylidene phosphorus ylide',
        explain: 'Correct — the C=C bond disconnects into an aldehyde and a phosphorus ylide (Ph–CH=PPh3), the classic Wittig retrosynthesis for alkenes.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(), C(), RING()]),
        label: 'Bibenzyl (alone)',
        explain: 'This is the saturated analogue — there is no simple substitution that installs a C=C bond from this starting material.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(), LBL('Br')]),
        mol2: MOL([RING(), C(), LBL('Br')]),
        condition: 'base',
        label: 'Two equivalents of benzyl bromide',
        explain: 'Simple deprotonation/alkylation of two benzyl halides does not form a C=C bond by any standard one-step method taught at this level.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C(), C(), LBL('H')], ['single', 'triple', 'single']),
        mol2: MOL([RING(), LBL('H')]),
        label: 'Phenylacetylene + Benzene',
        explain: 'This mixes an alkyne fragment with unfunctionalised benzene — there is no reasonable bond-forming step linking these to the target alkene.'
      }
    ]
  },
  {
    title: 'Acetophenone',
    target: MOL([RING(), C(BR('O', 'up', 'double')), C()]),
    highlightBond: 0,
    options: [
      {
        correct: true,
        mol1: MOL([RING()]),
        mol2: MOL([C(), C(BR('O', 'up', 'double')), LBL('Cl')]),
        condition: 'AlCl3',
        label: 'Benzene + Acetyl chloride',
        explain: 'Correct — disconnecting the aryl–carbonyl bond gives benzene and acetyl chloride: a Friedel–Crafts acylation.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('Cl')]),
        mol2: MOL([C(), C(BR('O', 'up', 'double')), LBL('OH')]),
        label: 'Chlorobenzene + Acetic acid',
        explain: 'Aryl halides do not undergo Friedel–Crafts acylation as the nucleophile in this way — the ring needs to be the nucleophile, not pre-halogenated.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C()]),
        label: 'Toluene (alone)',
        explain: 'This is the wrong oxidation pattern — oxidising the methyl group of toluene gives benzoic acid, not a ring-attached ketone.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('MgBr')]),
        mol2: MOL([C(), LBL('C≡N')]),
        label: 'Phenylmagnesium bromide + Acetonitrile',
        explain: 'This combination (Grignard + nitrile) is a more advanced route to the same product class, not the Friedel–Crafts disconnection being tested here.'
      }
    ]
  },
  {
    title: 'tert-Butylbenzene',
    target: MOL([RING(), LBL('C(CH3)3')]),
    highlightBond: 0,
    options: [
      {
        correct: true,
        mol1: MOL([RING()]),
        mol2: MOL([LBL('(CH3)3C'), LBL('Cl')]),
        condition: 'AlCl3',
        label: 'Benzene + tert-Butyl chloride',
        explain: 'Correct — disconnecting the aryl–alkyl bond gives benzene and tert-butyl chloride: a Friedel–Crafts alkylation, favoured because the tertiary carbocation is stable.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('Br')]),
        mol2: MOL([LBL('(CH3)3C'), LBL('MgBr')]),
        label: 'Bromobenzene + tert-Butylmagnesium bromide',
        explain: 'Two organometallic-style fragments do not couple directly under simple conditions — this is not the standard aromatic alkylation route.'
      },
      {
        correct: false,
        mol1: MOL([RING(), C()]),
        mol2: MOL([C(), C()]),
        label: 'Toluene + two methyl fragments',
        explain: 'There is no real bond-forming step here — "adding" separate methyl fragments to toluene is not a valid disconnection.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('OH')]),
        mol2: MOL([LBL('(CH3)3C'), LBL('Cl')]),
        label: 'Phenol + tert-Butyl chloride',
        explain: 'With phenol, the more nucleophilic oxygen would be alkylated instead, giving an ether — not the C-alkylated product shown.'
      }
    ]
  },
  {
    title: 'Anisole',
    target: MOL([RING(), LBL('O'), C()]),
    highlightBond: 1,
    options: [
      {
        correct: true,
        mol1: MOL([RING(), LBL('OH')]),
        mol2: MOL([C(), LBL('I')]),
        condition: 'K2CO3',
        label: 'Phenol + Methyl iodide',
        explain: 'Correct — disconnecting the O–CH3 bond gives phenol and methyl iodide: a Williamson ether synthesis (phenoxide attacks the methyl halide).'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('Cl')]),
        mol2: MOL([C(), LBL('OH')]),
        label: 'Chlorobenzene + Methanol',
        explain: 'Aryl chlorides are far too unreactive towards simple SN2 substitution by an alkoxide under standard conditions.'
      },
      {
        correct: false,
        mol1: MOL([RING()]),
        mol2: MOL([C(), LBL('OH')]),
        label: 'Benzene + Methanol',
        explain: 'This does not disconnect the marked O–CH3 bond at all, and there is no direct coupling method joining these two as written.'
      },
      {
        correct: false,
        mol1: MOL([RING(), LBL('OH')]),
        mol2: MOL([C(), LBL('OH')]),
        condition: 'H2SO4',
        label: 'Phenol + Methanol',
        explain: 'Simple acid-catalysed conditions do not effectively form aryl alkyl ethers this way — a good leaving group on the methyl partner is needed instead.'
      }
    ]
  }
];


/* ============================================================
   MAZE ENGINE — two linked rooms, progressive ghosts, live score
   Movement is discrete, cell-to-cell (see the maze prototype notes):
   an entity only ever interpolates between two pre-validated open
   cells, so leaving the maze is structurally impossible.
   ============================================================ */

var CELL = 32;
var DIRS = { up: { x: 0, y: -1 }, down: { x: 0, y: 1 }, left: { x: -1, y: 0 }, right: { x: 1, y: 0 } };
function opposite(d) { return { up: 'down', down: 'up', left: 'right', right: 'left' }[d]; }

var PLAYER_STEP_FRAMES = 20;
var GHOST_STEP_FRAMES = 27;

var GHOST_COLORS = ['#FF6B4A', '#00BFA6', '#B07CE8', '#FFD23F'];

/* ---- Room layouts ----
   '#' wall, '.' dot, '1'-'4' checkpoint (local index within the room,
   mapped to a global PUZZLES index via room.puzzleIndices), 'P' player
   start, 'G' a ghost present from the moment the room loads, 'X' the
   exit door (only usable once all 4 of the room's checkpoints are solved) */
var ROOMS = [
  {
    puzzleIndices: [0, 1, 2, 3], // Ethyl benzoate, Acetanilide, 1-Phenylethanol, Diphenylmethanol
    rows: [
      '#############',
      '#....G....G.#',
      '#.####.####.#',
      '#.#.1...2.#.#',
      '#.#.#####.#.#',
      '#.#.3...4.#.#',
      '#.####.####.#',
      '#.....P.....#',
      '######X######'
    ],
    // ghosts present from the start of this room
    ghostCount: 2,
    // additional ghosts that spawn partway through THIS room (none here —
    // the 3rd ghost is introduced on arrival in room 2 instead)
    bonusGhosts: []
  },
  {
    puzzleIndices: [4, 5, 6, 7], // (E)-Stilbene, Acetophenone, tert-Butylbenzene, Anisole
    rows: [
      '#############',
      '######X######',
      '#.....P.....#',
      '#.####.####.#',
      '#.#.1...2.#.#',
      '#.#.#####.#.#',
      '#.#.3...4.#.#',
      '#.####.####.#',
      '#.G........G#'
    ],
    ghostCount: 3, // 3rd ghost joins as you enter room 2 — difficulty ramps up
    // 4th ghost spawns once the 3rd checkpoint IN THIS ROOM is solved
    // (that's the 7th question overall)
    bonusGhosts: [{ afterLocalSolved: 3, col: 6, row: 3, color: GHOST_COLORS[3] }]
  }
];

var COLS = ROOMS[0].rows[0].length;
var ROWS = ROOMS[0].rows.length;

var grid, dots, checkpoints, doorCell, playerStart, ghostSpawns;

function parseRoom(roomIdx) {
  var room = ROOMS[roomIdx];
  grid = []; dots = []; checkpoints = {}; ghostSpawns = []; doorCell = null;
  for (var r = 0; r < room.rows.length; r++) {
    var rowArr = [], dotArr = [];
    for (var c = 0; c < room.rows[r].length; c++) {
      var ch = room.rows[r][c];
      if (ch === '#') { rowArr.push('#'); dotArr.push(false); continue; }
      if (ch === 'P') { playerStart = { col: c, row: r }; rowArr.push('.'); dotArr.push(true); continue; }
      if (ch === 'G') { ghostSpawns.push({ col: c, row: r }); rowArr.push('.'); dotArr.push(true); continue; }
      if (ch === 'X') { doorCell = { col: c, row: r }; rowArr.push('.'); dotArr.push(false); continue; }
      if (ch >= '1' && ch <= '4') {
        var localIdx = parseInt(ch, 10) - 1;
        checkpoints[r + ',' + c] = room.puzzleIndices[localIdx];
        rowArr.push('.'); dotArr.push(false); continue;
      }
      if (ch === '.') { rowArr.push('.'); dotArr.push(true); continue; }
      rowArr.push(' '); dotArr.push(false);
    }
    grid.push(rowArr); dots.push(dotArr);
  }
}

function isWallAt(col, row) {
  if (row < 0 || row >= grid.length || col < 0 || col >= grid[0].length) return true;
  return grid[row][col] === '#';
}
function cellCenter(col, row) { return { x: col * CELL + CELL / 2, y: row * CELL + CELL / 2 }; }

/* ---- Entities ---- */
var player, ghosts;

function makeGhost(col, row, color) {
  return { col: col, row: row, dir: null, moving: false, t: 0, color: color };
}

function resetPlayer() {
  player = { col: playerStart.col, row: playerStart.row, dir: null, desired: null, moving: false, t: 0, mouth: 0 };
}

function entityPixel(e) {
  var from = cellCenter(e.col, e.row);
  if (!e.moving) return from;
  var d = DIRS[e.dir];
  var to = cellCenter(e.col + d.x, e.row + d.y);
  return { x: from.x + (to.x - from.x) * e.t, y: from.y + (to.y - from.y) * e.t };
}

/* ---- Game / world state ---- */
var game = {
  room: 0,
  lives: 3,
  liveScore: 0,
  solved: [false, false, false, false, false, false, false, false],
  running: false,
  paused: false
};

function loadRoom(roomIdx) {
  game.room = roomIdx;
  parseRoom(roomIdx);
  resetPlayer();
  var room = ROOMS[roomIdx];
  ghosts = [];
  for (var i = 0; i < room.ghostCount; i++) {
    var s = ghostSpawns[i % ghostSpawns.length];
    ghosts.push(makeGhost(s.col, s.row, GHOST_COLORS[i]));
  }
  updateHud();
}

function localSolvedCountInRoom(roomIdx) {
  return ROOMS[roomIdx].puzzleIndices.filter(function (pi) { return game.solved[pi]; }).length;
}

function maybeSpawnBonusGhosts() {
  var room = ROOMS[game.room];
  room.bonusGhosts.forEach(function (bg) {
    var already = ghosts.some(function (g) { return g.color === bg.color; });
    if (!already && localSolvedCountInRoom(game.room) >= bg.afterLocalSolved) {
      ghosts.push(makeGhost(bg.col, bg.row, bg.color));
    }
  });
}

/* ---- Movement ---- */
function updatePlayer() {
  if (!player.moving) {
    if (dots[player.row][player.col]) {
      dots[player.row][player.col] = false;
      game.liveScore += 1;
      updateHud();
    }
    var key = player.row + ',' + player.col;
    if (checkpoints.hasOwnProperty(key) && !game.solved[checkpoints[key]]) {
      triggerCheckpoint(checkpoints[key], key);
      return;
    }
    if (doorCell && player.col === doorCell.col && player.row === doorCell.row) {
      if (localSolvedCountInRoom(game.room) === ROOMS[game.room].puzzleIndices.length) {
        advanceRoom();
        return;
      }
    }

    var order = [player.desired, player.dir].filter(Boolean);
    var chosen = null;
    for (var i = 0; i < order.length; i++) {
      var d = DIRS[order[i]];
      if (d && !isWallAt(player.col + d.x, player.row + d.y)) { chosen = order[i]; break; }
    }
    if (chosen) { player.dir = chosen; player.moving = true; player.t = 0; }
  } else {
    player.t += 1 / PLAYER_STEP_FRAMES;
    if (player.t >= 1) {
      var d2 = DIRS[player.dir];
      player.col += d2.x; player.row += d2.y;
      player.moving = false; player.t = 0;
    }
  }
  player.mouth = (player.mouth + 1) % 20;
}

function updateGhosts() {
  ghosts.forEach(function (g) {
    if (!g.moving) {
      var candidates = [];
      Object.keys(DIRS).forEach(function (name) {
        if (g.dir && name === opposite(g.dir)) return;
        var d = DIRS[name];
        if (!isWallAt(g.col + d.x, g.row + d.y)) candidates.push(name);
      });
      if (candidates.length === 0) {
        Object.keys(DIRS).forEach(function (name) {
          var d = DIRS[name];
          if (!isWallAt(g.col + d.x, g.row + d.y)) candidates.push(name);
        });
      }
      var best = candidates[0], bestDist = Infinity;
      candidates.forEach(function (name) {
        var d = DIRS[name];
        var nx = (g.col + d.x) - player.col, ny = (g.row + d.y) - player.row;
        var dist = nx * nx + ny * ny;
        if (dist < bestDist) { bestDist = dist; best = name; }
      });
      g.dir = best; g.moving = true; g.t = 0;
    } else {
      g.t += 1 / GHOST_STEP_FRAMES;
      if (g.t >= 1) {
        var d3 = DIRS[g.dir];
        g.col += d3.x; g.row += d3.y;
        g.moving = false; g.t = 0;
      }
    }
  });
}

function checkGhostCollision() {
  var pp = entityPixel(player);
  var hit = ghosts.some(function (g) {
    var gp = entityPixel(g);
    var dx = gp.x - pp.x, dy = gp.y - pp.y;
    return (dx * dx + dy * dy) < (CELL * 0.55) * (CELL * 0.55);
  });
  if (hit) loseLife();
}

function loseLife() {
  game.lives--;
  updateHud();
  if (game.lives <= 0) { endGame(false); return; }
  resetPlayer();
  ghostSpawns.forEach(function (s, i) {
    if (ghosts[i]) { ghosts[i].col = s.col; ghosts[i].row = s.row; ghosts[i].moving = false; ghosts[i].t = 0; ghosts[i].dir = null; }
  });
  // any bonus ghost beyond the base spawn points resets to its own defined spot
  ROOMS[game.room].bonusGhosts.forEach(function (bg) {
    var g = ghosts.filter(function (gg) { return gg.color === bg.color; })[0];
    if (g) { g.col = bg.col; g.row = bg.row; g.moving = false; g.t = 0; g.dir = null; }
  });
}

function advanceRoom() {
  if (game.room < ROOMS.length - 1) {
    loadRoom(game.room + 1);
  } else {
    endGame(true);
  }
}

/* ---- Checkpoint quiz overlay ---- */
var activeCheckpoint = null;

function triggerCheckpoint(pIdx, key) {
  game.paused = true;
  activeCheckpoint = { pIdx: pIdx, key: key };
  var puzzle = PUZZLES[pIdx];

  $('q-count').textContent = 'Question ' + (pIdx + 1) + ' of ' + PUZZLES.length;
  $('q-title').textContent = puzzle.title;
  $('q-target').innerHTML = renderMol(puzzle.target, { highlightBond: puzzle.highlightBond });

  var shuffled = shuffle(puzzle.options);
  var optsEl = $('q-options');
  optsEl.innerHTML = '';
  $('q-feedback').classList.add('hidden');

  shuffled.forEach(function (opt) {
    var btn = document.createElement('button');
    btn.className = 'option-card';
    btn.innerHTML = renderPrecursorPair(opt.mol1, opt.mol2) +
      (opt.condition ? '<span class="option-condition">' + opt.condition + '</span>' : '') +
      '<span class="option-condition"><strong>' + opt.label + '</strong></span>';
    btn.addEventListener('click', function () { answerCheckpoint(opt, btn, optsEl); });
    optsEl.appendChild(btn);
  });

  $('quiz-overlay').classList.remove('hidden');
  document.querySelector('.maze-wrap').classList.add('hidden');
  document.querySelector('.dpad').classList.add('hidden');
}

function answerCheckpoint(opt, btnEl, optsEl) {
  optsEl.querySelectorAll('.option-card').forEach(function (b) { b.disabled = true; });
  btnEl.classList.add(opt.correct ? 'correct' : 'incorrect');
  $('q-feedback-text').textContent = opt.explain;
  $('q-feedback').classList.remove('hidden');

  if (opt.correct) {
    game.solved[activeCheckpoint.pIdx] = true;
    delete checkpoints[activeCheckpoint.key];
    game.liveScore += 5;
  }
}

function closeCheckpointOverlay() {
  $('quiz-overlay').classList.add('hidden');
  document.querySelector('.maze-wrap').classList.remove('hidden');
  document.querySelector('.dpad').classList.remove('hidden');
  var solved = game.solved[activeCheckpoint.pIdx];
  activeCheckpoint = null;
  game.paused = false;
  maybeSpawnBonusGhosts();
  updateHud();
  if (!solved) {
    loseLife();
  }
}

/* ---- Rendering ---- */
var canvas, ctx;

// A small, subtle benzene-ring motif stamped on every wall tile — reuses
// the same hexagon-plus-circle convention as the molecule diagrams, at
// low opacity so it reads as texture rather than competing with the game.
function drawWallHex(cx, cy) {
  var r = 9;
  ctx.save();
  ctx.globalAlpha = 0.35;
  ctx.strokeStyle = '#7C8BFF';
  ctx.lineWidth = 1.1;
  ctx.beginPath();
  for (var k = 0; k < 6; k++) {
    var ang = (Math.PI / 180) * (60 * k - 90);
    var x = cx + r * Math.cos(ang), y = cy + r * Math.sin(ang);
    if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.closePath();
  ctx.stroke();
  ctx.beginPath();
  ctx.arc(cx, cy, r * 0.55, 0, Math.PI * 2);
  ctx.stroke();
  ctx.restore();
}

function drawMaze() {
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  for (var r = 0; r < grid.length; r++) {
    for (var c = 0; c < grid[0].length; c++) {
      var cx = c * CELL + CELL / 2, cy = r * CELL + CELL / 2;
      if (grid[r][c] === '#') {
        ctx.fillStyle = '#2B3A8F';
        ctx.fillRect(c * CELL + 2, r * CELL + 2, CELL - 4, CELL - 4);
        drawWallHex(cx, cy);
      } else if (dots[r][c]) {
        ctx.fillStyle = '#F5D76E';
        ctx.beginPath();
        ctx.arc(cx, cy, 2.6, 0, Math.PI * 2);
        ctx.fill();
      }
    }
  }

  // door (only visually "unlocked" once all of this room's checkpoints are solved)
  if (doorCell) {
    var unlocked = localSolvedCountInRoom(game.room) === ROOMS[game.room].puzzleIndices.length;
    var dcx = doorCell.col * CELL + CELL / 2, dcy = doorCell.row * CELL + CELL / 2;
    ctx.fillStyle = unlocked ? '#16A970' : '#4A4F6E';
    ctx.fillRect(dcx - 12, dcy - 12, 24, 24);
    ctx.fillStyle = '#10132A';
    ctx.font = 'bold 9px sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillText(unlocked ? 'GO' : 'LOCK', dcx, dcy + 1);
  }

  // checkpoints
  Object.keys(checkpoints).forEach(function (key) {
    var parts = key.split(','), r2 = +parts[0], c2 = +parts[1];
    var idx = checkpoints[key];
    var cx2 = c2 * CELL + CELL / 2, cy2 = r2 * CELL + CELL / 2;
    ctx.fillStyle = '#FF6B4A';
    ctx.beginPath(); ctx.arc(cx2, cy2, 10, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#10132A';
    ctx.font = 'bold 12px sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillText(String(idx + 1), cx2, cy2 + 1);
  });

  // ghosts
  ghosts.forEach(function (g) {
    var gp = entityPixel(g);
    ctx.fillStyle = g.color;
    ctx.beginPath();
    ctx.arc(gp.x, gp.y, 11, Math.PI, 0);
    ctx.lineTo(gp.x + 11, gp.y + 10);
    ctx.lineTo(gp.x + 5.5, gp.y + 5);
    ctx.lineTo(gp.x, gp.y + 10);
    ctx.lineTo(gp.x - 5.5, gp.y + 5);
    ctx.lineTo(gp.x - 11, gp.y + 10);
    ctx.closePath();
    ctx.fill();
    ctx.fillStyle = '#fff';
    ctx.beginPath(); ctx.arc(gp.x - 4, gp.y - 2, 3, 0, Math.PI * 2); ctx.fill();
    ctx.beginPath(); ctx.arc(gp.x + 4, gp.y - 2, 3, 0, Math.PI * 2); ctx.fill();
  });

  // player
  var pp = entityPixel(player);
  var open = (Math.abs(player.mouth - 10) / 10) * 0.28 + 0.04;
  var angleFor = { up: -Math.PI / 2, down: Math.PI / 2, left: Math.PI, right: 0 };
  var base = angleFor[player.dir || 'right'];
  ctx.fillStyle = '#F5D76E';
  ctx.beginPath();
  ctx.moveTo(pp.x, pp.y);
  ctx.arc(pp.x, pp.y, 11, base + open * Math.PI, base + (2 - open) * Math.PI);
  ctx.closePath();
  ctx.fill();
}

/* ---- HUD ---- */
function updateHud() {
  $('hud-lives').textContent = '♥ '.repeat(Math.max(game.lives, 0)).trim() || 'No lives left';
  var solvedCount = game.solved.filter(Boolean).length;
  $('hud-checkpoints').textContent = 'Solved ' + solvedCount + ' / ' + PUZZLES.length;
  $('hud-score').textContent = 'Score ' + game.liveScore;
}

/* ---- Main loop ---- */
function loop() {
  if (!game.running) return;
  if (!game.paused) {
    updatePlayer();
    if (!game.paused) {
      updateGhosts();
      checkGhostCollision();
    }
  }
  drawMaze();
  requestAnimationFrame(loop);
}

/* ---- Start / end screens ---- */
function shuffle(arr) {
  var a = arr.slice();
  for (var i = a.length - 1; i > 0; i--) {
    var j = Math.floor(Math.random() * (i + 1));
    var t = a[i]; a[i] = a[j]; a[j] = t;
  }
  return a;
}
function $(id) { return document.getElementById(id); }

function startGame() {
  game.lives = 3;
  game.liveScore = 0;
  game.solved = [false, false, false, false, false, false, false, false];
  game.paused = false;
  game.running = true;
  loadRoom(0);

  $('screen-intro').classList.add('hidden');
  $('screen-end').classList.add('hidden');
  $('screen-game').classList.remove('hidden');
  $('quiz-overlay').classList.add('hidden');
  document.querySelector('.maze-wrap').classList.remove('hidden');
  document.querySelector('.dpad').classList.remove('hidden');

  requestAnimationFrame(loop);
}

function endGame(won) {
  game.running = false;
  $('screen-game').classList.add('hidden');
  $('screen-end').classList.remove('hidden');

  var correctCount = game.solved.filter(Boolean).length;
  var pct = Math.round((correctCount / PUZZLES.length) * 100);
  var passed = won && game.lives > 0;

  $('end-heading').textContent = passed ? 'Maze cleared!' : 'Out of lives.';
  $('end-score').textContent = 'Checkpoints solved: ' + correctCount + ' / ' + PUZZLES.length +
    '. Lives remaining: ' + Math.max(game.lives, 0) + '. In-game score: ' + game.liveScore + ' (not sent to Moodle).';

  var reviewEl = $('end-review');
  reviewEl.innerHTML = '';
  PUZZLES.forEach(function (p, i) {
    var row = document.createElement('div');
    row.className = 'review-row ' + (game.solved[i] ? 'correct' : 'incorrect');
    row.innerHTML = '<span class="review-mark">' + (game.solved[i] ? '✓' : '✗') + '</span><span>' + p.title + '</span>';
    reviewEl.appendChild(row);
  });

  if (window.SCORM) {
    SCORM.setScore(pct, 100);
    SCORM.setComplete(passed);
  }
}

/* ---- Input ---- */
var KEY_MAP = {
  ArrowUp: 'up', ArrowDown: 'down', ArrowLeft: 'left', ArrowRight: 'right',
  w: 'up', s: 'down', a: 'left', d: 'right', W: 'up', S: 'down', A: 'left', D: 'right'
};

document.addEventListener('DOMContentLoaded', function () {
  canvas = $('maze');
  ctx = canvas.getContext('2d');

  if (window.SCORM) SCORM.init();

  $('btn-start').addEventListener('click', startGame);
  $('btn-restart').addEventListener('click', startGame);
  $('btn-q-continue').addEventListener('click', closeCheckpointOverlay);

  document.addEventListener('keydown', function (e) {
    var d = KEY_MAP[e.key];
    if (d) { player.desired = d; e.preventDefault(); }
  });

  [['dp-up', 'up'], ['dp-down', 'down'], ['dp-left', 'left'], ['dp-right', 'right']].forEach(function (pair) {
    var el = $(pair[0]);
    el.addEventListener('touchstart', function (e) { player.desired = pair[1]; e.preventDefault(); }, { passive: false });
    el.addEventListener('mousedown', function () { player.desired = pair[1]; });
  });

  window.addEventListener('beforeunload', function () {
    if (window.SCORM) SCORM.finish();
  });
});
