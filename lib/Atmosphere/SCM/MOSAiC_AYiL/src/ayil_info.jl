
"""
 The MOSAiC AYiL case set: the 190 days of the published 
 DALES dataset (Schnierstein et al., 2024, JAMES; Zenodo).
"""
const AYIL_DATES = (
    "20191016",
    "20191022",
    "20191024",
    "20191025",
    "20191026",
    "20191027",
    "20191028",
    "20191029",
    "20191030",
    "20191031",
    "20191101",
    "20191102",
    "20191103",
    "20191105",
    "20191106",
    "20191107",
    "20191108",
    "20191109",
    "20191110",
    "20191111",
    "20191112",
    "20191113",
    "20191114",
    "20191115",
    "20191125",
    "20191126",
    "20191129",
    "20191130",
    "20191201",
    "20191203",
    "20191204",
    "20191205",
    "20191206",
    "20191209",
    "20191210",
    "20191211",
    "20191212",
    "20191213",
    "20191214",
    "20191215",
    "20191217",
    "20191218",
    "20191219",
    "20191221",
    "20191222",
    "20191223",
    "20191225",
    "20191227",
    "20191229",
    "20191230",
    "20191231",
    "20200101",
    "20200102",
    "20200103",
    "20200105",
    "20200106",
    "20200109",
    "20200110",
    "20200111",
    "20200112",
    "20200114",
    "20200115",
    "20200116",
    "20200117",
    "20200118",
    "20200119",
    "20200120",
    "20200121",
    "20200122",
    "20200123",
    "20200124",
    "20200126",
    "20200127",
    "20200128",
    "20200130",
    "20200131",
    "20200202",
    "20200203",
    "20200204",
    "20200205",
    "20200206",
    "20200207",
    "20200208",
    "20200209",
    "20200210",
    "20200211",
    "20200212",
    "20200214",
    "20200215",
    "20200216",
    "20200217",
    "20200218",
    "20200219",
    "20200220",
    "20200221",
    "20200222",
    "20200223",
    "20200224",
    "20200225",
    "20200226",
    "20200227",
    "20200228",
    "20200301",
    "20200302",
    "20200303",
    "20200304",
    "20200305",
    "20200306",
    "20200307",
    "20200308",
    "20200310",
    "20200311",
    "20200316",
    "20200317",
    "20200318",
    "20200319",
    "20200321",
    "20200323",
    "20200324",
    "20200325",
    "20200326",
    "20200327",
    "20200329",
    "20200330",
    "20200402",
    "20200403",
    "20200404",
    "20200405",
    "20200406",
    "20200407",
    "20200408",
    "20200409",
    "20200410",
    "20200411",
    "20200412",
    "20200413",
    "20200414",
    "20200415",
    "20200416",
    "20200417",
    "20200418",
    "20200419",
    "20200420",
    "20200422",
    "20200423",
    "20200424",
    "20200425",
    "20200426",
    "20200427",
    "20200429",
    "20200430",
    "20200501",
    "20200502",
    "20200503",
    "20200504",
    "20200505",
    "20200506",
    "20200701",
    "20200702",
    "20200706",
    "20200707",
    "20200708",
    "20200709",
    "20200710",
    "20200713",
    "20200714",
    "20200715",
    "20200717",
    "20200720",
    "20200721",
    "20200722",
    "20200723",
    "20200724",
    "20200725",
    "20200726",
    "20200826",
    "20200827",
    "20200828",
    "20200829",
    "20200830",
    "20200831",
    "20200901",
    "20200902",
    "20200903",
    "20200905",
    "20200906",
    "20200907",
    "20200909",
    "20200910",
    "20200911",
)

# ============================================================================================= #




# ============================================================================================= #



"""
Best height [m] to simulate to to avoid anomalous ice.

Many MOSAiC AYiL days are archived with anomalous ice concentrations (much higher than reality,
as high as 10^8-10^8  per m³, due to initializing at 55 μm (described in the paper), which for SB3 
coefficients is much smaller than a typical size which is closer to 1mm.

These ice essentially never sediment or autoconvert, and behave very anomalously compared to ice
which develops naturally. Communication with Schnierstein indicates this may have been in error,
but initialization focused primarily on capturing the radiative impact of pre-existing ice, 
rather than on capturing the evolution of that ice.

However this makes it very challenging to guarantee that a reproduction of these simulations is
well-posed with different ice microphysicsal schemes and assumptions, ice-snow splits, etc.

Thus we provide a curated best-simulation top height for each day in this dictionary, for those cases
for which meaningful cloud exists beneath the lowest anomalous ice layer.

See (`docs/design.md` §12) for more details.

To guarantee easy reproduction, including ice outside the domain not causing sedimentation into domains
of interest not captured by the MOSAIC AYiL simulations, we recommend using these best heights.

The best split method we have found so far is a ice-size / fall-speed filer, provided in ___.
We apply to the full list, then apply a per-day hand check that
moved tops both up and down — 20200419's filter allows 4517 m but its cloud top
puts it on the 2500 m rung.

Finally, to avoid rampant recompilation and specialization, which massively bloats compilation times and RAM usage,
we condense onto three canonical grid heights, 2500m, 5000m, and LES_TOP_FACE.


These have been generated by <>, then hand checked



The value truncates the *domain*, not just the comparison: ice above the
trustworthy level would otherwise sediment and autoconvert down into the compared
region. 52 days run the full column; 19 are cut to 2500 or 5000 m.

A table rather than a computation, because the heights are not re-derivable. They
come from an ice-size and a fall-speed filter, t
"""

const RAW_BEST_SIMULATION_TOP_C = Dict{String, Float64}(
  "20200419" => 4517.39,
  "20200127" => 11857.2,
  "20200706" => 4823.92,
  "20200110" => 11857.2,
  "20200710" => 11857.2,
  "20200503" => 11857.2,
  "20200713" => 11857.2,
  "20191029" => 11857.2,
  "20200709" => 11857.2,
  "20200209" => 4166.48,
  "20200724" => 6177.97,
  "20200211" => 11857.2,
  "20200715" => 11857.2,
  "20191201" => 11857.2,
  "20200714" => 11857.2,
  "20191022" => 4248.64,
  "20191108" => 11857.2,
  "20200304" => 11857.2,
  "20200903" => 11857.2,
  "20191031" => 11857.2,
  "20191114" => 11857.2,
  "20191230" => 11857.2,
  "20200723" => 11857.2,
  "20200829" => 11857.2,
  "20191209" => 11857.2,
  "20200227" => 4823.92,
  "20200909" => 11857.2,
  "20200420" => 11857.2,
  "20200416" => 11857.2,
  "20200429" => 665.0,
  "20200911" => 11857.2,
  "20200418" => 11857.2,
  "20191218" => 11857.2,
  "20200831" => 11857.2,
  "20200826" => 11857.2,
  "20200103" => 11857.2,
  "20200827" => 2992.58,
  "20200707" => 3737.98,
  "20200409" => 3802.57,
  "20191219" => 4823.92,
  "20191213" => 6512.34,
  "20200905" => 11857.2,
  "20200124" => 3158.57,
  "20200828" => 11857.2,
  "20191221" => 11857.2,
  "20191225" => 11857.2,
  "20200726" => 11857.2,
  "20200502" => 11857.2,
  "20200901" => 3802.57,
  "20200224" => 11857.2,
  "20200210" => 2808.79,
  "20191210" => 11857.2,
  "20191103" => 11857.2,
  "20200307" => 11857.2,
  "20200410" => 4248.64,
  "20200708" => 11857.2,
  "20200907" => 11857.2,
  "20191125" => 11857.2,
  "20191016" => 11857.2,
  "20200830" => 11857.2,
  "20200425" => 11857.2,
  "20191028" => 6864.91,
  "20191030" => 11857.2,
  "20191024" => 11857.2,
  "20191101" => 11857.2,
  "20200216" => 3032.35,
  "20200415" => 4717.19,
  "20200717" => 875.0,
  "20191217" => 11857.2,
  "20200702" => 2273.43,
  "20200720" => 11857.2,
  "20200721" => 2099.65,
  "20200725" => 8157.23,
  "20200906" => 1983.08,
  "20191115" => 11857.2,
  "20200902" => 11857.2,
);

const RAW_BEST_SIMULATION_TOP_F  = Dict{String, Float64}(
  k => MOSAiC_AYiL_face_above_center(v) for (k,v) in RAW_BEST_SIMULATION_TOP_C)

const BEST_SIMULATION_TOP_F = Dict(date => RAW_BEST_SIMULATION_TOP_F[date] for date in keys(RAW_BEST_SIMULATION_TOP_F))


BEST_SIMULATION_TOP_F["20200210"] = 2500.0
BEST_SIMULATION_TOP_F["20200827"] = 2500.0
BEST_SIMULATION_TOP_F["20200216"] = 2500.0
BEST_SIMULATION_TOP_F["20200124"] = 2500.0
BEST_SIMULATION_TOP_F["20200707"] = 2500.0
BEST_SIMULATION_TOP_F["20200409"] = 2500.0
BEST_SIMULATION_TOP_F["20200901"] = 5000.0
BEST_SIMULATION_TOP_F["20200209"] = 5000.0
BEST_SIMULATION_TOP_F["20191022"] = 5000.0
BEST_SIMULATION_TOP_F["20200410"] = 5000.0
BEST_SIMULATION_TOP_F["20200419"] = 2500.0
BEST_SIMULATION_TOP_F["20200415"] = 5000.0
BEST_SIMULATION_TOP_F["20191219"] = 5000.0
BEST_SIMULATION_TOP_F["20200227"] = 5000.0
BEST_SIMULATION_TOP_F["20200706"] = 2500.0
BEST_SIMULATION_TOP_F["20200724"] = 2500.0
BEST_SIMULATION_TOP_F["20191213"] = 2500.0
BEST_SIMULATION_TOP_F["20191028"] = 2500.0
BEST_SIMULATION_TOP_F["20200725"] = 2500.0



"""Days a comparison against the reference is meaningful on, ascending."""
best_dates() = sort!(collect(keys(BEST_SIMULATION_TOP_F)))





