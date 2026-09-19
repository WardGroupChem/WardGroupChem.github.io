/* ============================================================
   MOLECULE RENDERER  (reused, unchanged, from the retrosynthesis
   disconnection quiz — see project technical-learnings)
   ============================================================ */

var BOND = 30;
var RING_R = 15;
var BRANCH_LEN = 24;
var FONT_SIZE = 13;
var DBL_OFFSET = 2.6;
var MIN_SEG = 7;

function RING() { return { t: 'ring' }; }
function C(branch) { return branch ? { t: 'c', branch: branch } : { t: 'c' }; }
function LBL(text) { return { t: 'label', text: text }; }
function BR(text, dir, bond) { return { text: text, dir: dir || 'up', bond: bond || 'single' }; }
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

function endPad(atom) {
  if (atom.t === 'ring') return RING_R;
  if (atom.t === 'label') return textHalfWidth(atom.text);
  return 0;
}

function computeSegSigns(atoms, bonds) {
  var n = atoms.length;
  var branchIdx = -1;
  for (var i = 0; i < n; i++) {
    if (atoms[i].t === 'c' && atoms[i].branch) { branchIdx = i; break; }
  }
  var s0 = 1;
  if (branchIdx > 0) {
    var required = (atoms[branchIdx].branch.dir === 'down') ? 1 : -1;
    var parityPow = ((branchIdx - 1) % 2 === 0) ? 1 : -1;
    s0 = required * parityPow;
  }
  var segSign = [];
  for (i = 0; i < n - 1; i++) segSign.push(s0 * ((i % 2 === 0) ? 1 : -1));

  bonds = bonds || [];
  for (i = 0; i < n - 1; i++) {
    if (bonds[i] === 'triple') {
      var s = (i - 1 >= 0) ? segSign[i - 1] : segSign[i];
      segSign[i] = s;
      if (i + 1 <= n - 2) {
        segSign[i + 1] = s;
        // Re-alternate every bond AFTER the forced linear run from scratch.
        // Only the two bonds immediately flanking the triple bond (the ones
        // attached to the sp carbons) are required to be collinear with it —
        // anything further down the chain must kink away normally. Without
        // this, a later bond can accidentally inherit the same sign as the
        // linear run purely by parity coincidence (same k%2), stretching it
        // into one unbroken over-length line instead of a proper zigzag.
        for (var j = i + 2; j < n - 1; j++) {
          segSign[j] = -segSign[j - 1];
        }
      }
    }
  }
  return segSign;
}

function drawBond(x1, y1, x2, y2, pad1, pad2, order, dashed) {
  order = order === true ? 'double' : (order === false ? 'single' : (order || 'single'));
  var dx = x2 - x1, dy = y2 - y1;
  var len = Math.sqrt(dx * dx + dy * dy) || 1;
  var ux = dx / len, uy = dy / len;
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

  var pos = [{ x: 0, y: 0 }];
  for (var i = 0; i < n - 1; i++) {
    var pad1 = endPad(atoms[i]), pad2 = endPad(atoms[i + 1]);
    var segLen = Math.max(BOND, pad1 + pad2 + MIN_SEG + 4);
    var sdx, sdy;
    if (segSign[i] === 0) {
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

  for (i = 0; i < n - 1; i++) {
    var order = bonds[i] || 'single';
    var isHi = (i === highlightBond);
    body += drawBond(pos[i].x, pos[i].y, pos[i + 1].x, pos[i + 1].y, endPad(atoms[i]), endPad(atoms[i + 1]), order, isHi);
  }

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

/* Forward-synthesis reaction row: start molecule -> [reagents go here] -> target molecule.
   (Replaces the disconnection quiz's renderPrecursorPair, since this quiz asks
   "what reagents complete this synthesis" rather than "what are the two precursors". */
function renderReaction(startMol, targetMol) {
  var html = '<div class="reaction-row">';
  html += '<div class="reaction-mol"><span class="mol-cap">Start</span>' + renderMol(startMol) + '</div>';
  html += '<div class="reaction-arrow"><span class="arrow-q">?</span><span class="arrow-line"><span class="arrow-shaft"></span><span class="arrow-head"></span></span></div>';
  html += '<div class="reaction-mol"><span class="mol-cap">Target</span>' + renderMol(targetMol) + '</div>';
  html += '</div>';
  return html;
}

/* ============================================================
   PUZZLE DATA — alkene & alkyne synthesis, 2nd year undergrad
   ============================================================ */

var PUZZLES = [
  {
    title: 'But-2-ene, from an alkyl halide',
    start: MOL([C(), C(BR('Br', 'up')), C(), C()]),
    target: MOL([C(), C(), C(), C()], ['single', 'double', 'single']),
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
        explain: 'Br2 adds across a C=C double bond — it\'s a test/trap for an alkene, not a reagent that installs one from an alkyl halide.'
      }
    ]
  },
  {
    title: 'But-2-ene, from an alcohol',
    start: MOL([C(), C(BR('OH', 'up')), C(), C()]),
    target: MOL([C(), C(), C(), C()], ['single', 'double', 'single']),
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
    title: '(E)-Stilbene, by Wittig olefination',
    start: MOL([RING(), C(BR('O', 'up', 'double'))]),
    target: MOL([RING(), C(), C(), RING()], ['single', 'double', 'single']),
    options: [
      {
        correct: true,
        reagents: 'PPh₃, PhCH₂Br, then n-BuLi',
        explain: 'Correct — PPh3 and benzyl bromide form a phosphonium salt; deprotonation with base generates the benzylidene ylide, which reacts with the aldehyde shown to form the new C=C bond directly (a Wittig reaction).'
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
    start: MOL([C(), C(), C(), C(), C()], ['single', 'triple', 'single', 'single']),
    // Each branch sits exactly 120° from the double bond's own direction at
    // that carbon (same rule as the carbonyl branches) — both branches on
    // the same 120°-derived side for a cis (Z) alkene.
    target: MOL([C3([BR3('CH3', -90)]), C3([BR3('CH2CH3', -30)])], ['double']),
    options: [
      {
        correct: true,
        reagents: 'H₂, Lindlar catalyst',
        explain: 'Correct — Lindlar\'s poisoned Pd catalyst delivers H2 syn across the triple bond in a single step and stops at the alkene, giving the cis (Z) product.'
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
    start: MOL([C(), C(), C(), C(), C()], ['single', 'triple', 'single', 'single']),
    // Same 120°-from-the-double-bond rule, but the second branch takes the
    // OTHER 120°-valid position (opposite side) — trans (E) alkene.
    target: MOL([C3([BR3('CH3', -90)]), C3([BR3('CH2CH3', 90)])], ['double']),
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
    start: MOL([C(BR('Br', 'up')), C(BR('Br', 'down')), C(), C()]),
    target: MOL([C(), C(), C(), C()], ['triple', 'single', 'single']),
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
    start: MOL([C3([BR3('CH3', -135), BR3('CH3', -45), BR3('Br', 45), BR3('Br', 135)])]),
    target: MOL([C(), C(), C()], ['triple', 'single']),
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
    start: MOL([C(), C(), C(), C(), C()], ['triple', 'single', 'single', 'single']),
    target: MOL([C(), C(), C(), C(), C(), C(), C()], ['single', 'single', 'triple', 'single', 'single', 'single']),
    options: [
      {
        correct: true,
        reagents: 'LDA, THF; then CH₃CH₂Br',
        explain: 'Correct — LDA deprotonates the terminal alkyne to give the acetylide, a strong carbanion nucleophile, which performs SN2 substitution on the primary ethyl bromide to extend the chain.'
      },
      {
        correct: false,
        reagents: 'LDA, THF; then (CH₃)₃CBr',
        explain: 'The acetylide is a strong, hindered base — with a tertiary halide like tert-butyl bromide it promotes E2 elimination instead of SN2 substitution, so the chain isn\'t extended as drawn.'
      },
      {
        correct: false,
        reagents: 'CH₃CH₂Br (no base first)',
        explain: 'Without first forming the acetylide, there\'s no nucleophile present to attack the alkyl bromide — nothing happens.'
      },
      {
        correct: false,
        reagents: 'NaOH; then CH₃CH₂Br',
        explain: 'Hydroxide (conjugate acid pKa ~15.7) is nowhere near basic enough to deprotonate a terminal alkyne (pKa ~25) — no acetylide forms, so no alkylation can occur.'
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
