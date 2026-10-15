# Minnesota Glacial Till Analysis Skill(s) for Claude

This repo (currently) holds one skill intended for use by the Minnesota Geological Survey (MGS) in assisting their analysis
of glacial tills using texture and grain estimates.

## Skill: Review existing members and formations
The `review-till-formations` skill looks to analyze samples already coded by the geologists at MGS for statistical and geological
patterns and inconsistencies.  It does that through the following workflow

**1) Identify Bad Data:**
- Locate and remove samples that have errors such as null values, incorrect texture values (don't sum to 100), 
bad depth readings (negative or unreasonably large), and more
- Partially or fully remove these samples from analysis depending on what useful information can be extracted from them

**2) Create a geological timeline from existing data:**
- Come to a consensus as to **how** formations and members should be ordered chronologically (downhole ordering)
- Look for evidence through formation contacts in deep holes
- Boreholes can **skip** over a formation, but never reorder them
- Work with the user to solve ambiguous or heavily contradicting orderings

**3) Identify timeline violations:**
- Identify and filter coded samples that do **not** follow the timeline defined in step 2

**4) Identify sample outliers:**
- Identify outliers (<1% chance of occurring) 
- Can be caused by low overall grain count, a geological process (see `geological_priors`), or human error
- Attempt to give a cause for each outlier, and remove some or all of its data from analysis depending on how much useful 
information can be extracted from it

**5) Compute current formation distributions:**
- Compute the sand/silt/clay, precambrian/paleozoic/cretaceous, crystalline/carbonate/shale, and light/dark/red composition
of each formation within the county and of neighboring counties

**6) Test samples and contacts within each borehole:**
- Identify samples that are **most likely** to be from a different formation than they were coded as
- Test where the most probable contact range is between different formations in boreholes, and how it differs from the
existing contact
- Test different rates of formation mixture in potentially mixed samples to see if it better explains an ambiguous sample
as compared to a single code

**7) Test adjacent units to see if they are a single unit**
- Check if chronologically adjacent units are more likely to have come from a single distribution rather than two
- Report how strong the evidence is in terms of contacts and boreholes containing the formations

**8) Test units to see if they can be split into two (or more) formations**
- Check to see if a unit is more likely to be split into an upper and lower section than a single distribution
- Report how many boreholes both the upper and lower sections exist in, and how much vertical evidence there is versus lateral spread

### Statistical and Geological Priors

This skill looks to incorporate a mixture of geological priors to identify outlier caused by a variety of processes along with
dirichlet multinomial distributions to describe grain counts and bivariate normal to describe texture.

**Geological Priors:**
- Helps determine if odd looking data can be explained by leeching, sediment incorporation, weathering, etc.
- Encodes commonly misidentified or mislabeled grains and what they are confused for
- Enforces a timeline that must be followed

**Statistical methods:**
- Uses a dirichlet multinomial distribution to explain how samples are a draw from an unknown distribution that shifts vertically
  (over time) within a borehole and laterally (over space).
- Texture can be explained by a simple bivariate normal distribution