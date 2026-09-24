/* ============================================================
   Disconnect// maze — game logic.
   Structures: chem.js (SMILES -> SVG). Puzzle data: puzzles.js
   (shared with the standalone quiz, so both always match).
   ============================================================ */

/* ============================================================
   MAZE ENGINE — two linked rooms, progressive ghosts, live score
   Movement is discrete, cell-to-cell (see the maze prototype notes):
   an entity only ever interpolates between two pre-validated open
   cells, so leaving the maze is structurally impossible.
   ============================================================ */

var CELL = 32;
var DIRS = { up: { x: 0, y: -1 }, down: { x: 0, y: 1 }, left: { x: -1, y: 0 }, right: { x: 1, y: 0 } };
function opposite(d) { return { up: 'down', down: 'up', left: 'right', right: 'left' }[d]; }

var BASE_PLAYER_STEP_FRAMES = 24;
var BASE_GHOST_STEP_FRAMES = 32;
var playerStepFrames = BASE_PLAYER_STEP_FRAMES;
var ghostStepFrames = BASE_GHOST_STEP_FRAMES;

// Ghosts occasionally break off their chase for a short "scatter" burst of
// semi-random movement, so they feel less mechanically perfect.
/* ---- Difficulty presets (chosen on the intro screen) ----
   Player speed is the same on every setting so the controls always feel
   the same; difficulty only changes the ghosts.
     ghostSpeed  multiplier on ghost frames-per-cell (>1 = slower ghosts)
     ghostDelta  ghosts removed from each room's starting count
     bonusGhost  whether the extra ghost joins after question 7
     flee        multiplier on how long catalyst pellets make ghosts flee
     scatter     chance a ghost breaks into random wandering at a junction */
var DIFFICULTY_PRESETS = {
  easy:   { label: 'Easy',   ghostSpeed: 1.25, ghostDelta: 1, bonusGhost: false, flee: 1.5,  scatter: 0.22 },
  normal: { label: 'Normal', ghostSpeed: 1.0,  ghostDelta: 0, bonusGhost: true,  flee: 1.0,  scatter: 0.15 },
  hard:   { label: 'Hard',   ghostSpeed: 0.85, ghostDelta: 0, bonusGhost: true,  flee: 0.75, scatter: 0.08 }
};
var difficulty = 'normal';
function preset() { return DIFFICULTY_PRESETS[difficulty] || DIFFICULTY_PRESETS.normal; }

var SCATTER_CHANCE = 0.15; // overwritten from the preset at game start
var SCATTER_DURATION_FRAMES = BASE_GHOST_STEP_FRAMES * 3;

// Catalyst pellets make every ghost flee (and become eatable) for this long.
var FLEE_DURATION_FRAMES = 480; // ~8s at ~60fps
var FLEE_FLASH_THRESHOLD_FRAMES = 90; // start flashing in the last ~1.5s

var GHOST_COLORS = ['#FF6B4A', '#00BFA6', '#B07CE8', '#FFD23F'];

/* ---- Room layouts ----
   '#' wall, '.' dot, '1'-'4' checkpoint (local index within the room,
   mapped to a global PUZZLES index via a per-attempt shuffled order —
   see loadRoom), 'P' player start, 'G' a ghost present from the moment
   the room loads, 'X' the exit door (only usable once all 4 of the
   room's checkpoints are solved), 'C' a catalyst power-up pellet. */
var ROOMS = [
  {
    puzzleIndices: [0, 1, 2, 3], // Ethyl benzoate, Acetanilide, 1-Phenylethanol, Diphenylmethanol
    rows: [
      '#############',
      '#.C..G....G.#',
      '#.####.####.#',
      '#.#.1...2.#.#',
      '#.#.#####.#.#',
      '#.#.3...4.#.#',
      '#.####.####.#',
      '#.....P...C.#',
      '######X######'
    ],
    // ghosts present from the start of this room
    ghostCount: 2,
    // additional ghosts that spawn partway through THIS room (none here —
    // the 3rd ghost is introduced on arrival in room 2 instead)
    bonusGhosts: []
  },
  {
    // A deliberately different topology from room 1: a three-rail "ladder"
    // (verticals at columns 2, 6 and 10 threading through four horizontal
    // corridors) instead of room 1's twin side-boxes, so a player who
    // remembers room 1's shape doesn't get a free pass here.
    puzzleIndices: [4, 5, 6, 7], // (E)-Stilbene, Acetophenone, tert-Butylbenzene, Anisole
    rows: [
      '#############',
      '#G.1.....2.G#',
      '##.###.###.##',
      '#...........#',
      '##.###.###.##',
      '#.G.3.C.4...#',
      '##.###.###.##',
      '#.....P..C..#',
      '######X######'
    ],
    ghostCount: 3, // 3rd ghost joins as you enter room 2 — difficulty ramps up
    // 4th ghost spawns once the 3rd checkpoint IN THIS ROOM is solved
    // (that's the 7th question overall)
    bonusGhosts: [{ afterLocalSolved: 3, col: 6, row: 1, color: GHOST_COLORS[3] }]
  }
];

var COLS = ROOMS[0].rows[0].length;
var ROWS = ROOMS[0].rows.length;

var grid, dots, pellets, checkpoints, checkpointPositions, doorCell, playerStart, ghostSpawns;

// roomIdx: which room to parse. order: the shuffled array of global PUZZLES
// indices for this room's 4 checkpoints (order[0] -> checkpoint '1', etc) —
// re-rolled by loadRoom() on every attempt so replay layouts aren't fixed.
function parseRoom(roomIdx, order) {
  var room = ROOMS[roomIdx];
  grid = []; dots = []; pellets = {}; checkpoints = {}; checkpointPositions = {};
  ghostSpawns = []; doorCell = null;
  for (var r = 0; r < room.rows.length; r++) {
    var rowArr = [], dotArr = [];
    for (var c = 0; c < room.rows[r].length; c++) {
      var ch = room.rows[r][c];
      if (ch === '#') { rowArr.push('#'); dotArr.push(false); continue; }
      if (ch === 'P') { playerStart = { col: c, row: r }; rowArr.push('.'); dotArr.push(true); continue; }
      if (ch === 'G') { ghostSpawns.push({ col: c, row: r }); rowArr.push('.'); dotArr.push(true); continue; }
      if (ch === 'X') { doorCell = { col: c, row: r }; rowArr.push('.'); dotArr.push(false); continue; }
      if (ch === 'C') { pellets[r + ',' + c] = true; rowArr.push('.'); dotArr.push(false); continue; }
      if (ch >= '1' && ch <= '4') {
        var localIdx = parseInt(ch, 10) - 1;
        var pIdx = order[localIdx];
        checkpoints[r + ',' + c] = pIdx;
        checkpointPositions[r + ',' + c] = pIdx;
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
  // homeCol/homeRow are fixed at creation and never move — every ghost,
  // however it was spawned (base or bonus), resets to exactly this point
  // when a life is lost. This is what the old build got wrong: it only
  // knew how to reset the first couple of ghosts.
  return { col: col, row: row, homeCol: col, homeRow: row, dir: null, moving: false, t: 0, color: color, scatterUntil: 0 };
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
  paused: false,
  muted: true,
  frameCount: 0,
  fleeUntilFrame: 0
};

function applyRoomSpeed(roomIdx) {
  if (roomIdx === 0) {
    playerStepFrames = BASE_PLAYER_STEP_FRAMES;
    ghostStepFrames = Math.round(BASE_GHOST_STEP_FRAMES * preset().ghostSpeed);
  } else {
    // ~10% faster (fewer frames per cell) — a small escalation on top of
    // the extra ghost count, so room 2 reads as a genuine step up.
    playerStepFrames = Math.round(BASE_PLAYER_STEP_FRAMES * 0.9);
    ghostStepFrames = Math.round(BASE_GHOST_STEP_FRAMES * 0.9 * preset().ghostSpeed);
  }
}

function loadRoom(roomIdx) {
  game.room = roomIdx;
  var room = ROOMS[roomIdx];
  var order = shuffle(room.puzzleIndices.slice());
  parseRoom(roomIdx, order);
  resetPlayer();
  ghosts = [];
  // spawns are listed farthest-first in each room, so Easy drops the
  // ghost that starts closest to the player
  var nGhosts = Math.max(1, room.ghostCount - preset().ghostDelta);
  for (var i = 0; i < nGhosts; i++) {
    var s = ghostSpawns[i % ghostSpawns.length];
    ghosts.push(makeGhost(s.col, s.row, GHOST_COLORS[i]));
  }
  applyRoomSpeed(roomIdx);
  game.fleeUntilFrame = 0;
  updateHud();
}

function localSolvedCountInRoom(roomIdx) {
  return ROOMS[roomIdx].puzzleIndices.filter(function (pi) { return game.solved[pi]; }).length;
}

function maybeSpawnBonusGhosts() {
  var room = ROOMS[game.room];
  if (!preset().bonusGhost) return;
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
      playSound('chomp');
      updateHud();
    }
    var key = player.row + ',' + player.col;
    if (pellets[key]) {
      delete pellets[key];
      game.liveScore += 3;
      game.fleeUntilFrame = game.frameCount + Math.round(FLEE_DURATION_FRAMES * preset().flee);
      playSound('powerUp');
      updateHud();
    }
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
    player.t += 1 / playerStepFrames;
    if (player.t >= 1) {
      var d2 = DIRS[player.dir];
      player.col += d2.x; player.row += d2.y;
      player.moving = false; player.t = 0;
    }
  }
  player.mouth = (player.mouth + 1) % 20;
}

function updateGhosts() {
  var fleeing = game.frameCount < game.fleeUntilFrame;
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

      var best;
      if (fleeing) {
        // run AWAY from the player — maximise distance instead of minimising it
        var worstDist = -Infinity;
        candidates.forEach(function (name) {
          var d = DIRS[name];
          var nx = (g.col + d.x) - player.col, ny = (g.row + d.y) - player.row;
          var dist = nx * nx + ny * ny;
          if (dist > worstDist) { worstDist = dist; best = name; }
        });
      } else if (game.frameCount < (g.scatterUntil || 0)) {
        // mid-scatter: wander semi-randomly instead of chasing
        best = candidates[Math.floor(Math.random() * candidates.length)];
      } else {
        var bestDist = Infinity;
        candidates.forEach(function (name) {
          var d = DIRS[name];
          var nx = (g.col + d.x) - player.col, ny = (g.row + d.y) - player.row;
          var dist = nx * nx + ny * ny;
          if (dist < bestDist) { bestDist = dist; best = name; }
        });
        // small chance to break into a scatter burst at the next junction
        if (Math.random() < SCATTER_CHANCE) g.scatterUntil = game.frameCount + SCATTER_DURATION_FRAMES;
      }
      g.dir = best; g.moving = true; g.t = 0;
    } else {
      g.t += 1 / ghostStepFrames;
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
  var fleeing = game.frameCount < game.fleeUntilFrame;
  for (var i = 0; i < ghosts.length; i++) {
    var g = ghosts[i];
    var gp = entityPixel(g);
    var dx = gp.x - pp.x, dy = gp.y - pp.y;
    if ((dx * dx + dy * dy) < (CELL * 0.55) * (CELL * 0.55)) {
      if (fleeing) {
        g.col = g.homeCol; g.row = g.homeRow;
        g.moving = false; g.t = 0; g.dir = null;
        game.liveScore += 10;
        updateHud();
        playSound('eatGhost');
      } else {
        playSound('ghostCatch');
        loseLife();
        return;
      }
    }
  }
}

function loseLife() {
  game.lives--;
  updateHud();
  if (game.lives <= 0) { endGame(false); return; }
  resetPlayer();
  // every ghost (base spawn or bonus) returns to the exact point it was
  // created at — this is the fix for the old "3rd ghost never resets" bug
  ghosts.forEach(function (g) {
    g.col = g.homeCol; g.row = g.homeRow;
    g.moving = false; g.t = 0; g.dir = null;
    g.scatterUntil = 0;
  });
  game.fleeUntilFrame = 0;
}

function advanceRoom() {
  if (game.room < ROOMS.length - 1) {
    playSound('roomClear');
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
  $('q-target').innerHTML = Chem.svg(puzzle.target, { scale: 1.2, title: puzzle.title });

  var shuffled = shuffle(puzzle.options);
  var optsEl = $('q-options');
  optsEl.innerHTML = '';
  $('q-feedback').classList.add('hidden');

  shuffled.forEach(function (opt) {
    var btn = document.createElement('button');
    btn.className = 'option-card';
    btn.innerHTML = renderPrecursors(opt) +
      (opt.condition ? '<span class="option-condition">' + chemText(opt.condition) + '</span>' : '') +
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
  $('q-feedback-text').innerHTML = chemText(opt.explain);
  $('q-feedback').classList.remove('hidden');

  if (opt.correct) {
    game.solved[activeCheckpoint.pIdx] = true;
    delete checkpoints[activeCheckpoint.key];
    game.liveScore += 5;
    playSound('correct');
  } else {
    playSound('incorrect');
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

/* ---- Sound (synthesised, no audio files; muted by default) ---- */
var audioCtx = null;
function getAudioCtx() {
  if (!audioCtx) {
    var AC = window.AudioContext || window.webkitAudioContext;
    if (!AC) return null;
    audioCtx = new AC();
  }
  return audioCtx;
}

function beep(freq, durationMs, type, startDelay, gainLevel) {
  var ctx2 = getAudioCtx();
  if (!ctx2) return;
  var t0 = ctx2.currentTime + (startDelay || 0);
  var osc = ctx2.createOscillator();
  var gain = ctx2.createGain();
  osc.type = type || 'square';
  osc.frequency.setValueAtTime(freq, t0);
  gain.gain.setValueAtTime(0, t0);
  gain.gain.linearRampToValueAtTime(gainLevel || 0.15, t0 + 0.01);
  gain.gain.linearRampToValueAtTime(0, t0 + durationMs / 1000);
  osc.connect(gain); gain.connect(ctx2.destination);
  osc.start(t0);
  osc.stop(t0 + durationMs / 1000 + 0.02);
}

function playSound(name) {
  if (game.muted) return;
  switch (name) {
    case 'chomp': beep(660, 40, 'square', 0, 0.05); break;
    case 'correct': beep(523, 90, 'sine', 0, 0.15); beep(784, 140, 'sine', 0.09, 0.15); break;
    case 'incorrect': beep(220, 180, 'sawtooth', 0, 0.15); break;
    case 'eatGhost': beep(880, 60, 'square', 0, 0.15); beep(1174, 90, 'square', 0.05, 0.15); break;
    case 'powerUp': beep(392, 60, 'triangle', 0, 0.15); beep(523, 60, 'triangle', 0.06, 0.15); beep(659, 90, 'triangle', 0.12, 0.15); break;
    case 'ghostCatch': beep(160, 220, 'sawtooth', 0, 0.18); break;
    case 'roomClear': beep(523, 100, 'sine', 0, 0.15); beep(659, 100, 'sine', 0.1, 0.15); beep(784, 100, 'sine', 0.2, 0.15); beep(1047, 180, 'sine', 0.3, 0.15); break;
    case 'gameOver': beep(392, 150, 'sawtooth', 0, 0.15); beep(311, 150, 'sawtooth', 0.15, 0.15); beep(220, 300, 'sawtooth', 0.3, 0.15); break;
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

  // catalyst pellets — larger than a dot and gently pulsing so they read
  // as a distinct, deliberate pickup rather than just a bigger dot
  var pulse = 4 + Math.sin(game.frameCount * 0.15) * 1.4;
  Object.keys(pellets).forEach(function (key) {
    var parts = key.split(','), r3 = +parts[0], c3 = +parts[1];
    var cx3 = c3 * CELL + CELL / 2, cy3 = r3 * CELL + CELL / 2;
    ctx.fillStyle = '#00E5C7';
    ctx.beginPath();
    ctx.arc(cx3, cy3, pulse, 0, Math.PI * 2);
    ctx.fill();
  });

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

  // checkpoints — unsolved ones show as numbered orange markers; solved
  // ones leave a faint green tick behind instead of vanishing outright,
  // so the maze itself shows your progress at a glance
  Object.keys(checkpointPositions).forEach(function (key) {
    var parts = key.split(','), r2 = +parts[0], c2 = +parts[1];
    var idx = checkpointPositions[key];
    var cx2 = c2 * CELL + CELL / 2, cy2 = r2 * CELL + CELL / 2;
    if (checkpoints.hasOwnProperty(key)) {
      ctx.fillStyle = '#FF6B4A';
      ctx.beginPath(); ctx.arc(cx2, cy2, 10, 0, Math.PI * 2); ctx.fill();
      ctx.fillStyle = '#10132A';
      ctx.font = 'bold 12px sans-serif';
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(String(idx + 1), cx2, cy2 + 1);
    } else {
      ctx.save();
      ctx.globalAlpha = 0.55;
      ctx.strokeStyle = '#16A970';
      ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(cx2, cy2, 8, 0, Math.PI * 2); ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(cx2 - 4, cy2);
      ctx.lineTo(cx2 - 1, cy2 + 3);
      ctx.lineTo(cx2 + 4, cy2 - 4);
      ctx.stroke();
      ctx.restore();
    }
  });

  // ghosts — pale blue (flashing white near the end) while fleeing/eatable
  var fleeing = game.frameCount < game.fleeUntilFrame;
  var flashOn = fleeing && (game.fleeUntilFrame - game.frameCount) < FLEE_FLASH_THRESHOLD_FRAMES &&
    (Math.floor(game.frameCount / 8) % 2 === 0);
  ghosts.forEach(function (g) {
    var gp = entityPixel(g);
    ctx.fillStyle = fleeing ? (flashOn ? '#FFFFFF' : '#5AC8FA') : g.color;
    ctx.beginPath();
    ctx.arc(gp.x, gp.y, 11, Math.PI, 0);
    ctx.lineTo(gp.x + 11, gp.y + 10);
    ctx.lineTo(gp.x + 5.5, gp.y + 5);
    ctx.lineTo(gp.x, gp.y + 10);
    ctx.lineTo(gp.x - 5.5, gp.y + 5);
    ctx.lineTo(gp.x - 11, gp.y + 10);
    ctx.closePath();
    ctx.fill();
    ctx.fillStyle = fleeing ? '#10132A' : '#fff';
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
  $('hud-room').textContent = 'Room ' + (game.room + 1) + '/' + ROOMS.length + ' · ' + preset().label;
  $('hud-lives').textContent = '♥ '.repeat(Math.max(game.lives, 0)).trim() || 'No lives left';
  var solvedCount = game.solved.filter(Boolean).length;
  $('hud-checkpoints').textContent = 'Solved ' + solvedCount + ' / ' + PUZZLES.length;
  $('hud-score').textContent = 'Score ' + game.liveScore;
}

/* ---- Main loop ---- */
function loop() {
  if (!game.running) return;
  if (!game.paused) {
    game.frameCount++;
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
  game.frameCount = 0;
  game.fleeUntilFrame = 0;
  SCATTER_CHANCE = preset().scatter;
  loadRoom(0);

  // AudioContext needs a user gesture to start — this click qualifies
  var ac = getAudioCtx();
  if (ac && ac.state === 'suspended') ac.resume();

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
  playSound(passed ? 'roomClear' : 'gameOver');

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
  // "Try again" returns to the intro so the difficulty can be changed
  $('btn-restart').addEventListener('click', function () {
    $('screen-end').classList.add('hidden');
    $('screen-intro').classList.remove('hidden');
    window.scrollTo(0, 0);
  });

  document.querySelectorAll('.diff-btn').forEach(function (btn) {
    btn.addEventListener('click', function () {
      difficulty = btn.getAttribute('data-diff');
      document.querySelectorAll('.diff-btn').forEach(function (b) {
        b.classList.toggle('active', b === btn);
        b.setAttribute('aria-pressed', b === btn ? 'true' : 'false');
      });
    });
  });
  $('btn-q-continue').addEventListener('click', closeCheckpointOverlay);

  var muteBtn = $('hud-mute');
  muteBtn.textContent = game.muted ? '🔇' : '🔊';
  muteBtn.addEventListener('click', function () {
    game.muted = !game.muted;
    muteBtn.textContent = game.muted ? '🔇' : '🔊';
    if (!game.muted) {
      var ac2 = getAudioCtx();
      if (ac2 && ac2.state === 'suspended') ac2.resume();
    }
  });

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
