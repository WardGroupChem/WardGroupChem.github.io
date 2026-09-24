/* ============================================================
   PUZZLE DATA — shared by the quiz and the maze.
   Every structure is a SMILES string drawn by chem.js, so bond
   lengths/angles/labels are consistent by construction.
   Local SMILES extensions (see chem.js):
     ^         before a bond  -> that bond is the dashed red disconnection
     {A|B}     abbreviation label (A when bonded on its left, B when
               bonded on its right), digits drawn as subscripts
   ============================================================ */

var PUZZLES = [
  {
    title: "Ethyl benzoate",
    target: "CCO^C(=O)c1ccccc1",
    options: [
      {
        correct: true,
        mol1: "OC(=O)c1ccccc1",
        mol2: "CCO",
        condition: "H⁺, heat",
        label: "Benzoic acid + Ethanol",
        explain: "Correct — disconnecting the acyl–oxygen bond gives benzoic acid and ethanol. Forward direction: Fischer esterification (acid + alcohol, H⁺ catalyst)."
      },
      {
        mol1: "OCc1ccccc1",
        mol2: "CC(=O)O",
        label: "Benzyl alcohol + Acetic acid",
        explain: "The acyl and alkyl portions are swapped — this pair would give benzyl acetate, not ethyl benzoate."
      },
      {
        mol1: "Oc1ccccc1",
        mol2: "CCC(=O)O",
        label: "Phenol + Propanoic acid",
        explain: "Wrong connectivity and chain length — this would give phenyl propanoate, a different ester entirely."
      },
      {
        mol1: "ClC(=O)c1ccccc1",
        mol2: "CCCO",
        label: "Benzoyl chloride + Propan-1-ol",
        explain: "The acyl chloride is a reasonable acylating agent, but the alcohol has one carbon too many — this gives propyl benzoate."
      }
    ]
  },
  {
    title: "Acetanilide",
    target: "CC(=O)^Nc1ccccc1",
    options: [
      {
        correct: true,
        mol1: "Nc1ccccc1",
        mol2: "CC(=O)Cl",
        condition: "base",
        label: "Aniline + Acetyl chloride",
        explain: "Correct — disconnecting the N–C(=O) bond gives aniline and acetyl chloride. Forward direction: acylation of the amine."
      },
      {
        mol1: "NC(=O)c1ccccc1",
        mol2: "CN",
        label: "Benzamide + Methylamine",
        explain: "This breaks a different bond and swaps which fragment carries the carbonyl — it leads to N-methylbenzamide, an isomeric but different amide."
      },
      {
        mol1: "{NO2|O2N}c1ccccc1",
        mol2: "CC(=O)Cl",
        label: "Nitrobenzene + Acetyl chloride",
        explain: "Nitrobenzene has the wrong oxidation state at nitrogen — it would need reduction to aniline first before this disconnection applies."
      },
      {
        mol1: "Oc1ccccc1",
        mol2: "CC(N)=O",
        label: "Phenol + Acetamide",
        explain: "This pair of functional groups does not combine to form the marked N–C bond — it points toward an ester, not an amide."
      }
    ]
  },
  {
    title: "1-Phenylethanol",
    target: "C^C(O)c1ccccc1",
    options: [
      {
        correct: true,
        mol1: "O=Cc1ccccc1",
        mol2: "C{MgBr|BrMg}",
        condition: "then H₃O⁺",
        label: "Benzaldehyde + Methylmagnesium bromide",
        explain: "Correct — the C–CH3 bond next to the alcohol is disconnected as a Grignard addition: the methyl nucleophile adds to the aldehyde carbonyl."
      },
      {
        mol1: "CC(=O)c1ccccc1",
        label: "Acetophenone (alone)",
        explain: "This is the oxidised ketone, not a disconnection — reducing it would work synthetically, but it does not break the marked C–C bond retrosynthetically."
      },
      {
        mol1: "{MgBr|BrMg}c1ccccc1",
        mol2: "CC=O",
        label: "Phenylmagnesium bromide + Acetaldehyde",
        explain: "This breaks the Ph–C bond instead of the marked C–CH3 bond — a valid Grignard route to the same product, but not the disconnection shown here."
      },
      {
        mol1: "CCc1ccccc1",
        label: "Ethylbenzene (alone)",
        explain: "There is no bond-forming logic here — this is simply the deoxygenated hydrocarbon, not a retrosynthetic precursor pair."
      }
    ]
  },
  {
    title: "Diphenylmethanol",
    target: "OC(c1ccccc1)^c1ccccc1",
    options: [
      {
        correct: true,
        mol1: "O=Cc1ccccc1",
        mol2: "{MgBr|BrMg}c1ccccc1",
        condition: "then H₃O⁺",
        label: "Benzaldehyde + Phenylmagnesium bromide",
        explain: "Correct — disconnecting the marked C–Ph bond gives benzaldehyde and a phenyl Grignard, which adds to the carbonyl."
      },
      {
        mol1: "O=C(c1ccccc1)c1ccccc1",
        label: "Benzophenone (alone)",
        explain: "This is the oxidised ketone. Reduction (e.g. NaBH4) gets you to the product, but it is not a disconnection of the marked bond."
      },
      {
        mol1: "c1ccccc1Cc1ccccc1",
        label: "Diphenylmethane (alone)",
        explain: "This lacks the oxygen entirely and has no bond-forming step that reconnects to give the alcohol."
      },
      {
        mol1: "c1ccc(cc1)-c1ccccc1",
        mol2: "[H]C([H])=O",
        label: "Biphenyl + Formaldehyde",
        explain: "Wrong connectivity — biphenyl already has the two rings joined directly, which is not the bond pattern in the target."
      }
    ]
  },
  {
    title: "(E)-Stilbene",
    target: "c1ccccc1/C^=C/c1ccccc1",
    options: [
      {
        correct: true,
        mol1: "O=Cc1ccccc1",
        mol2: "c1ccccc1C={PPh3|Ph3P}",
        condition: "Wittig",
        label: "Benzaldehyde + a benzylidene phosphorus ylide",
        explain: "Correct — the C=C bond disconnects into an aldehyde and a phosphorus ylide (Ph–CH=PPh3), the classic Wittig retrosynthesis for alkenes."
      },
      {
        mol1: "c1ccccc1CCc1ccccc1",
        label: "Bibenzyl (alone)",
        explain: "This is the saturated analogue — there is no simple substitution that installs a C=C bond from this starting material."
      },
      {
        mol1: "BrCc1ccccc1",
        mol2: "BrCc1ccccc1",
        condition: "base",
        label: "Two equivalents of benzyl bromide",
        explain: "Simple deprotonation/alkylation of two benzyl halides does not form a C=C bond by any standard one-step method taught at this level."
      },
      {
        mol1: "[H]C#Cc1ccccc1",
        mol2: "c1ccccc1",
        label: "Phenylacetylene + Benzene",
        explain: "This mixes an alkyne fragment with unfunctionalised benzene — there is no reasonable bond-forming step linking these to the target alkene."
      }
    ]
  },
  {
    title: "Acetophenone",
    target: "CC(=O)^c1ccccc1",
    options: [
      {
        correct: true,
        mol1: "c1ccccc1",
        mol2: "CC(=O)Cl",
        condition: "AlCl3",
        label: "Benzene + Acetyl chloride",
        explain: "Correct — disconnecting the aryl–carbonyl bond gives benzene and acetyl chloride: a Friedel–Crafts acylation."
      },
      {
        mol1: "Clc1ccccc1",
        mol2: "CC(=O)O",
        label: "Chlorobenzene + Acetic acid",
        explain: "Aryl halides do not undergo Friedel–Crafts acylation as the nucleophile in this way — the ring needs to be the nucleophile, not pre-halogenated."
      },
      {
        mol1: "Cc1ccccc1",
        label: "Toluene (alone)",
        explain: "This is the wrong oxidation pattern — oxidising the methyl group of toluene gives benzoic acid, not a ring-attached ketone."
      },
      {
        mol1: "{MgBr|BrMg}c1ccccc1",
        mol2: "CC#N",
        label: "Phenylmagnesium bromide + Acetonitrile",
        explain: "This combination (Grignard + nitrile) is a more advanced route to the same product class, not the Friedel–Crafts disconnection being tested here."
      }
    ]
  },
  {
    title: "tert-Butylbenzene",
    target: "CC(C)(C)^c1ccccc1",
    options: [
      {
        correct: true,
        mol1: "c1ccccc1",
        mol2: "CC(C)(C)Cl",
        condition: "AlCl3",
        label: "Benzene + tert-Butyl chloride",
        explain: "Correct — disconnecting the aryl–alkyl bond gives benzene and tert-butyl chloride: a Friedel–Crafts alkylation, favoured because the tertiary carbocation is stable."
      },
      {
        mol1: "Brc1ccccc1",
        mol2: "CC(C)(C){MgBr|BrMg}",
        label: "Bromobenzene + tert-Butylmagnesium bromide",
        explain: "Two organometallic-style fragments do not couple directly under simple conditions — this is not the standard aromatic alkylation route."
      },
      {
        mol1: "Cc1ccccc1",
        mol2: "[CH3]",
        condition: "× 2",
        label: "Toluene + two methyl fragments",
        explain: "There is no real bond-forming step here — \"adding\" separate methyl fragments to toluene is not a valid disconnection."
      },
      {
        mol1: "Oc1ccccc1",
        mol2: "CC(C)(C)Cl",
        label: "Phenol + tert-Butyl chloride",
        explain: "With phenol, the more nucleophilic oxygen would be alkylated instead, giving an ether — not the C-alkylated product shown."
      }
    ]
  },
  {
    title: "Anisole",
    target: "C^Oc1ccccc1",
    options: [
      {
        correct: true,
        mol1: "Oc1ccccc1",
        mol2: "CI",
        condition: "K2CO3",
        label: "Phenol + Methyl iodide",
        explain: "Correct — disconnecting the O–CH3 bond gives phenol and methyl iodide: a Williamson ether synthesis (phenoxide attacks the methyl halide)."
      },
      {
        mol1: "Clc1ccccc1",
        mol2: "CO",
        label: "Chlorobenzene + Methanol",
        explain: "Aryl chlorides are far too unreactive towards simple SN2 substitution by an alkoxide under standard conditions."
      },
      {
        mol1: "c1ccccc1",
        mol2: "CO",
        label: "Benzene + Methanol",
        explain: "This does not disconnect the marked O–CH3 bond at all, and there is no direct coupling method joining these two as written."
      },
      {
        mol1: "Oc1ccccc1",
        mol2: "CO",
        condition: "H2SO4",
        label: "Phenol + Methanol",
        explain: "Simple acid-catalysed conditions do not effectively form aryl alkyl ethers this way — a good leaving group on the methyl partner is needed instead."
      }
    ]
  }
];

/* Formula text (explanations, conditions): digits after an element symbol
   or a closing bracket become subscripts, e.g. CH3, NaBH4, K2CO3. */
function chemText(str) {
  var esc = String(str).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  return esc.replace(/([A-Za-z\)])(\d+)/g, '$1<sub>$2</sub>');
}

/* One option's starting materials, drawn at a common scale and height. */
function renderPrecursors(opt) {
  var html = '<div class="option-precursors">';
  if (opt.mol2) {
    var pair = Chem.pair(opt.mol1, opt.mol2);
    html += pair[0] + '<span class="plus-sign">+</span>' + pair[1];
  } else {
    html += Chem.svg(opt.mol1);
  }
  return html + '</div>';
}
