# Geological priors reference

Every geological prior decision in a report must come from a section below. Some of these priors have methods associated 
with them implemented in `scripts/geological_priors/`.

Contents: §1 Outliers · §2 Grain count human error

## §1 Outliers 

A sample can sit outside its unit expected distribution because a **process** changed its grain counts or texture, not because it belongs to a 
different unit. Each process below lists what it alters, the columns that identify it, and what separates it from a 
real change of unit.

**Rule common to all of them.** A process acts on some parts of the composition and leaves the rest alone, and it is 
tied to a position (the surface, a contact, bedrock, sorted sediment). So a process is accepted only when:
1. the altered parts moved in the expected direction against the sample's neighbors in the same hole,
2. the composition **with those parts removed** is typical of those neighbors, and
3. an independent indicator (flag, colour, depth, adjacent deposit, note) is present.

1 + 2 without 3 is "possible". 3 alone means the process was present and changed nothing measurable.

**Indicator columns**

| Column                                                       | Use | Caveat                                                |
|--------------------------------------------------------------|---|-------------------------------------------------------|
| `leached`                                                    | Y = geologist saw no HCl reaction | blank on a large % of till samples; blank is not N    |
| `wet_color`, `dry_color`                                     | Munsell chroma ≥ 3 = oxidized, ≤ 2 = reduced | dry color much less commomn than wet                  |
| `sampledepth`, `depth_top`, `depth_bottom`                   | distance from surface; interval sampled | `depth_unit = m` must be converted to feet first      |
| `deposittype` of rows above/below in the same hole           | O outwash, L lacustrine, B bedrock next to the sample | blank on some rows                                    |
| `secondary`, `unknown`, `other_other`                        | authigenic and non-provenance grains | what they are is only in the notes                    |
| `pyrite`, `lignite`                                          | unstable in oxidized till |                                                       |
| `sand_12_pct`, `gravel_pct`, `silt_norm`                     | coarse tail and silt share of fines |                                                       |
| `notes1`–`5`, `grain_note1`–`4`, `previnterpt`, `strat_work` | free text, keywords listed per process | long notes are cut across columns: search them joined |

### §1.A Weathering (after deposition, from a surface downward)

**§1.1 Leaching.** Carbonate dissolves from the top of a till downward.
- *Grains:* `carbonate`, `limestone`, `inoceramus` fall, to zero when complete. Calcite (limestone, inoceramus) goes 
  before dolomite, so a partly leached sample can keep `carbonate` and have lost the two Cretaceous types. Every other 
  type rises in % because the total closes to 100; their counts relative to each other do not change.
- *Texture:* silt and clay fall slightly (fine carbonate is lost), sand rises.
- *Identify:* `leached = Y`; carbonate types at or near 0 while the samples below in the hole have them; carbonate 
  rising with depth across consecutive samples; notes `leach`, `noncalc`, `no carbonate`. Usually oxidized as well.
- *Not a different unit if:* the carbonate-free composition matches the neighbors. Zero carbonate alone is not 
  evidence: some tills (Superior, Rainy provenance) have almost none.
- *Buried leaching:* the top of an older till under a younger one can be leached (a weathering surface). Same 
  signature, but at depth and directly below a contact. It marks the top of a unit.

**§1.2 Oxidation.** Oxygenated water reaches the till above the water table and along joints and sand bodies.
- *Grains:* `pyrite` is destroyed and becomes iron oxide; `gray_shale`, `speckled_shale` and `lignite` soften and break 
  down, so they fall; `secondary` rises (limonite, goethite, Fe-cemented grains).
- *Texture:* little; clay may rise where shale has broken down.
- *Identify:* wet chroma ≥ 3 (hues 10YR, 7.5YR); `pyrite = 0` where reduced samples below have it; `secondary` above 
  the neighbors; notes `oxid`, `secondary Fe`, `limonite`, `goethite`, `iron oxide`, `stain`, `coat`.
- *Not a different unit if:* the composition without shale, pyrite and lignite matches the neighbors. Oxidized color 
  with typical counts is an indicator only.
- *At depth:* oxidized color far below the surface points to a buried weathering surface, a sand body next to the 
  sample, a joint, or oxidation of the core after drilling. Check `deposittype` of adjacent rows and the notes.

**§1.3 Shale breakdown (slaking).** Soft shale falls apart with wetting, drying and freezing, and again in the lab 
during washing and sieving.
- *Grains:* shale types fall with no change in the others. *Texture:* clay and silt rise by the same material.
- *Identify:* shale low and clay high against the neighbors together; notes `crumbl`, `broken up`, `disintegrat`, 
  `soft`. Shallow or oxidized position supports it. Lab breakdown has no positional indicator: only the notes.

**§1.4 Secondary minerals.** New grains grow in the till: secondary carbonate (caliche, concretions, root casts, 
calcite cement), siderite, Fe and Mn oxides, gypsum, pyrite in reduced till.
- *Grains:* `secondary` rises. Provenance counts are unchanged, but fewer provenance grains are counted. If secondary 
  carbonate is counted as `carbonate`, carbonate is too high; cemented till fragments land in `other_other`.
- *Texture:* cemented aggregates behave as sand: sand and `sand_12_pct` rise, clay falls.
- *Identify:* `secondary` as a share of all grains (above 10 % in 1,161 till samples, above 25 % in 276); `pyrite`; 
  notes `secondary`, `siderite`, `concretion`, `caliche`, `cement`, `gypsum`, `selenite`, `mangan`, `till ball`, 
  `till frag`.
- *Not a different unit if:* provenance counts alone match the neighbors.

**§1.5 Soil formation.** In the top few feet, clay moves from the A to the B horizon, roots and animals mix the 
profile, and organic matter is added. Usually with §1.1 and §1.2.
- *Texture:* clay low in the top sample and high in the one below, against deeper samples of the same unit.
- *Grains:* wood, roots and debris in `other_other`; otherwise as §1.1–1.2.
- *Identify:* `sampledepth` within about 5 ft of the surface; notes `soil`, `horizon`, `Bt`, `root`, `organic`, 
  `wood`. A buried soil (`paleosol`) gives the same at depth, below a contact.

**§1.6 Loess or other silt mixed into the top.** Windblown silt is worked into the upper till.
- *Texture:* silt and `silt_norm` rise, sand and gravel fall. *Grains:* unchanged, since loess has no 1–2 mm grains, 
  but `sand_12_pct` and the grain total are low.
- *Identify:* top sample of the hole; silt high with typical grains; notes `loess`, `silt cap`.

**§1.7 Weathering of crystalline grains.** In old or long-exposed tills, feldspar-bearing and mafic grains rot.
- *Grains:* `felsic`, `mafic_igneous`, `metasedvol` fall; the resistant types `clear_quartz`, `quartzite`, `chert` 
  rise as a share of Precambrian grains. *Texture:* clay rises.
- *Identify:* resistant share of Precambrian above the neighbors, together with §1.1 and §1.2; notes `weather`, 
  `rotten`, `punky`, `grus`, `saprol`, `kaolin`.
- *Caution:* felsic against clear quartz is also where grain counters disagree most (§2). Check `counted_by` first.

### §1.B Glacial and depositional processes (same ice, different local input)

**§1.8 Bedrock incorporation.** The base of a till picks up the rock or saprolite under it.
- *Grains:* the local lithology rises: Cretaceous (`gray_shale`, `speckled_shale`, `lignite`, `ostrandersand`, 
  `limestone`, `pyrite`), Paleozoic (`carbonate`, `paleo_ss`, `paleo_shale`), or weathered crystalline rock. 
  *Texture:* toward the bedrock: clay over shale, sand over sandstone.
- *Identify:* a `deposittype = B` row or a bedrock `stratunit` code within 10 ft below; the enrichment is largest in 
  the lowest sample and fades upward; notes `bedrock`, `saprol`, `local`.
- *Not a different unit if:* the composition without the local lithology matches the samples above.

**§1.9 Incorporation of older till.** Ice mixes the till under it into the base of the new one.
- *Grains and texture:* between the unit above and the unit below, in every part at once.
- *Identify:* the sample is at a contact, and its composition lies on the mixing line between the two units (a 
  weighted mean of the two in proportions, not in log-ratios); notes `mixed`, `incorporat`, `rework`.
- *Not a different unit if:* one mixing fraction fits grains and texture together. A sample that fits neither unit and 
  no mixture of them is a candidate thin unit.

**§1.10 Incorporation of sorted sediment.** Ice overrides outwash or lake sediment, or a sample includes a sand lens.
- *Texture:* over outwash, sand and gravel rise; over lake sediment, silt and clay rise and gravel and `sand_12_pct` 
  fall. *Grains:* typical of the unit, or slightly toward the overridden sediment; the grain total is low in the 
  fine case.
- *Identify:* `deposittype = O` or `L` within 5 ft in the same hole; texture shifted and grains typical; notes `sand 
  lens`, `stringer`, `laminat`, `interbed`, `sandy till`.

**§1.11 Washing and surface reworking.** Meltwater, waves (glacial lake plains), slope movement or frost remove 
fines from the top of a till or move it.
- *Texture:* clay and silt fall, sand and gravel rise (a gravel lag). *Grains:* typical.
- *Identify:* top of a unit or of the hole, or directly under sorted sediment; notes `washed`, `lag`, `colluv`, 
  `slopewash`, `flow till`, `frost`, `cryo`.

### §1.C Not geological: rule these out first

They produce the same numbers as the processes above and are checked from the record alone.
- *Small count:* fewer than 100 provenance grains; notes `not enough`, `too few`, `less than 200`.
- *Non-provenance grains:* `unknown` + `secondary` + `other_other` above 25 % of grains; notes `junk`, `debris`.
- *Sample spans a contact:* `depth_top`/`depth_bottom` interval crosses a change of `deposittype` or unit. Looks like 
  §1.9.
- *Drilling contamination:* material from higher in the hole in a deeper sample; notes `slough`, `cave`, 
  `contamin`, `cuttings`. Looks like §1.9 or §1.12.
- *Counter and recount:* `counted_by` differs from the neighbors; notes `recount`. See §2.
- *Record errors:* texture not summing to 100 ± 2, depth missing or in metres, duplicate depth.

## §2 Grain count human error

A sample can sit outside its unit expected distribution due to human error caused by similar looking grains.  Below is a 
list of the most commonly confused grains.

- `felsic` + `clear_quartz`
- `felsic` + `carbonate`
- `mafic_igneous` + `metasedvol`
- `inoceramus` + any
- `red_volcanic` + any