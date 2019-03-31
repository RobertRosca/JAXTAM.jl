# JAXTAM

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://RobertRosca.github.io/JAXTAM.jl/latest) [![Build Status](https://travis-ci.org/RobertRosca/JAXTAM.jl.svg?branch=master)](https://travis-ci.org/RobertRosca/JAXTAM.jl) [![Coverage Status](https://coveralls.io/repos/RobertRosca/JAXTAM.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/RobertRosca/JAXTAM.jl?branch=master) [![codecov.io](http://codecov.io/github/RobertRosca/JAXTAM.jl/coverage.svg?branch=master)](http://codecov.io/github/RobertRosca/JAXTAM.jl?branch=master)

Just Another X-ray Timing Analysis Module

JAXTAM, inspired by [MaLTPyNT](https://maltpynt.readthedocs.io/en/latest/), is a Julia package for X-ray timing analysis, with a specific focus on HEASARC ran missions.

## Install

1. Download Julia version 1.1.0 [from the website](https://julialang.org/downloads/), extract the archive to a convenient place and add the julia executable to your path\*

2. Clone this repository

3. Start julia

4. Press `]` to enter package mode, then type in `activate path_to_JAXTAM`

5. Write `instantiate` (still in package mode)

6. Press backspace to exit package mode

7. Type `using JAXTAM` to import the module, this will precompile the dependencies and may take a while

NOTE: if you are running this headless, on a server to generate reports, then you need to enable headless plotting, to do this
edit `~/.julia/config/startup.jl` and add in:

```
@info "Running headless plots:  'ENV[\"GKSwstype\"] = \"100\"'"
ENV["GKSwstype"] = "100"
```

## Basic Usage

### 1. Path Setup

First, you need to set up the required paths:

* download - Directory data will be downloaded to
* jaxtam - Directory JAXTAM-created data files will be saved to
* web - Directory plots and HTML report pages are saved to
* rmf - Path to the mission's RMF file

To do this run `JAXTAM.mission_paths(mission_name)`, for NICER this looks like:

```
julia> JAXTAM.mission_paths(nicer)
[ Info: Mission not found in /home/robertr/Projects/JAXTAM.jl/mission_paths.json, please enter paths:
Download path: /export/data/robertr/heasarc/nicer/download
JAXTAM (processed data) path: /export/data/robertr/heasarc/nicer/jaxtam  
Web (html reports) path: /export/data/robertr/heasarc/nicer/web     
RMF (caldb mission file) path: /home/sw-astro/caldb/data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf
┌ Info: Wrote to /home/robertr/Projects/JAXTAM.jl/mission_paths.json
└       Add custom keys in to JSON file manually if required
(download = "/export/data/robertr/heasarc/nicer/download", jaxtam = "/export/data/robertr/heasarc/nicer/jaxtam", web = "/export/data/robertr/heasarc/nicer/web", rmf = "/home/sw-astro/caldb/data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf")
```

### 2. Master Table Setup

Now you need to download the master table for the mission, call `JAXTAM.master(mission)`:

```
julia> JAXTAM.master(nicer)
[ Info: Downloading latest master catalog
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  664k  100  664k    0     0   330k      0  0:00:02  0:00:02 --:--:--  330k

7-Zip [64] 9.20  Copyright (c) 1999-2010 Igor Pavlov  2010-11-18
p7zip Version 9.20 (locale=en_GB.UTF-8,Utf16=on,HugeFiles=on,40 CPUs)

Processing archive: /export/data/robertr/heasarc/nicer/jaxtam/master.tdat

Extracting  master

Everything is Ok

Size:       4955750
Compressed: 680884
[ Info: Loading /export/data/robertr/heasarc/nicer/jaxtam/master.tdat
[ Info: Saving /export/data/robertr/heasarc/nicer/jaxtam/master.feather
[ Info: Loading /export/data/robertr/heasarc/nicer/jaxtam/master.feather
[ Info: Generating append table, this may take some time
[ Info: Checking publicity
[ Info: Looping through logs
        Log 1248/1248
```

Depending on how many observations have been analysed/downloaded the append table generation can take a few minutes

### 3. Observation Querying

Now that the tables are set up, you can use query functions to pick which observation(s) you want to analyse

The basic syntax for this is `JAXTAM.master_query(mission, :key, value)`, for example looking for all observations by name:

```
julia> JAXTAM.master_query(nicer, :name, "X_Persei")
4×35 DataFrames.DataFrame. Omitted printing of 28 columns
│ Row │ name     │ ra      │ dec     │ lii     │ bii      │ time                │ end_time            │
│     │ String   │ Float64 │ Float64 │ Float64 │ Float64  │ Dates.DateTime      │ Dates.DateTime      │
├─────┼──────────┼─────────┼─────────┼─────────┼──────────┼─────────────────────┼─────────────────────┤
│ 1   │ X_Persei │ 58.851  │ 31.0455 │ 163.085 │ -17.1337 │ 2019-01-01T21:45:40 │ 2019-01-01T22:26:20 │
│ 2   │ X_Persei │ 58.8458 │ 31.0457 │ 163.081 │ -17.1365 │ 2018-12-10T01:00:30 │ 2018-12-10T01:41:00 │
│ 3   │ X_Persei │ 58.8458 │ 31.0457 │ 163.081 │ -17.1365 │ 2018-12-08T16:35:11 │ 2018-12-08T17:16:40 │
│ 4   │ X_Persei │ 58.8454 │ 31.0459 │ 163.081 │ -17.1366 │ 2018-12-09T00:18:10 │ 2018-12-09T13:20:20 │
```

If you're interested in only publicly available observations, then use `master_query_public` instead

Additionally, just calling `master_query_public(mission)` will return all of the public observations:

```
julia> size(JAXTAM.master(nicer), 1)
11172

julia> size(JAXTAM.master_query_public(nicer), 1)
10986
```

Above we see there are 11,172 observations with 10986 of those being marked as public

The more practical search would be for public, non-calibration, observations, which requires some filtering:

```
julia> size(filter(o->o[:obs_type]!="CAL", JAXTAM.master_query_public(nicer)), 1)
9867
```

So we see there are 9,867 public non-calibration observations

### 4. Observation Downloading

The download commands take in either `DataFrame` rows (as returned by the query functions) or `obsid`s

For example, if you know you want to download `1200360101` all you need to do is:

```
julia> JAXTAM.download(nicer, "1200360101")
[ Info: heasarc.gsfc.nasa.gov:/.nicer_archive/.nicer_201809a/obs/2018_09/1200360101 --> /export/data/robertr/heasarc/nicer/download/nicer_archive/nicer_201809a/obs/2018_09/1200360101
`lftp heasarc.gsfc.nasa.gov -e 'mirror "/.nicer_archive/.nicer_201809a/obs/2018_09/1200360101" "/export/data/robertr/heasarc/nicer/download/nicer_archive/nicer_201809a/obs/2018_09/1200360101" --parallel=10 --only-newer --exclude-glob *ufa.evt.gz --exclude-glob *ufa.evt --exclude-glob *uf.evt.gz && exit'`
Total: 6 directories, 14 files, 0 symlinks
```

Or if you want to download all public non-calibration observations:

```
julia> download_queue = filter(o->o[:obs_type]!="CAL", JAXTAM.master_query_public(nicer))
11226×35 DataFrames.DataFrame. Omitted printing of 29 columns, 11226 rows
julia> JAXTAM.download(nicer, download_queue)
...
```

### 5. Observation Reports

As with the downloading, report generation can take in an `obsid` or a table of observations

If you want to see help for a function, type in `?` to enter help mode, then the name of the function:

```
help?> JAXTAM.report
  report(mission::Mission, obs_row::DataFrameRow; e_range::Tuple{Float64,Float64}=_mission_good_e_range(mission), overwrite::Bool, nuke::Bool=false, update_masterpage::Bool=true)

  Generates a report for the default energy range

  Creates plots for:

  * Lightcurve (+ grouped lightcurves)
  
  * Periodigram (+ grouped periodograms)
  
  * Power Spectra
  
      * :rms, full range, log-rebinned, log-log plot
  
      * :leahy, full range, log-rebinned, log-log plot
  
      * :leahy, 0 to 1 Hz, no rebin, linear-linear plot
  
      * :leahy, 1 to end Hz, no rebin, linear-linear plot
  
      * :leahy, 50 to end Hz, no regin, linear-linear plot
  
  * Spectrogram
  
  * Pulsation search plot

  Produces HTML report page

  Updates the homepage

  ──────────────────────────────────────────────────────────────────────────────────────────────────────────

  report(mission::Mission, obsid::String; kwargs...)

  Multiple dispath to report(mission::Mission, obs_row::DataFrameRow; kwargs...)

  Uses obsid to select the observation

```

As you can see, there are two ways to use this function: either give it `obs_rows` from a query result, or an `obsid`

If you want to create reports for multiple observations, then basic list comprehension is the easiest way to do this:

```
julia> report_queue = filter(o->o[:obs_type]!="CAL", JAXTAM.master_query_public(nicer, :name, "YZ_CMi"))
3×35 DataFrames.DataFrame. Omitted printing of 28 columns
│ Row │ name   │ ra      │ dec     │ lii     │ bii     │ time                │ end_time            │
│     │ String │ Float64 │ Float64 │ Float64 │ Float64 │ Dates.DateTime      │ Dates.DateTime      │
├─────┼────────┼─────────┼─────────┼─────────┼─────────┼─────────────────────┼─────────────────────┤
│ 1   │ YZ_CMi │ 116.17  │ 3.55356 │ 215.856 │ 13.4601 │ 2019-01-28T04:14:50 │ 2019-01-28T08:03:00 │
│ 2   │ YZ_CMi │ 116.169 │ 3.55426 │ 215.855 │ 13.4598 │ 2019-01-27T03:53:26 │ 2019-01-27T08:52:00 │
│ 3   │ YZ_CMi │ 116.169 │ 3.55428 │ 215.855 │ 13.4596 │ 2019-01-26T04:20:51 │ 2019-01-26T08:08:20 │
julia> [JAXTAM.report(nicer, obs) for obs in eachrow(report_queue)]
```

By default reports are generated for the full mission energy range as returned by `JAXTAM._mission_good_e_range(mission)`

You can set this to an alternate energy range with the keyword `e_range` as a tuple of floats (in keV): `e_range=(0.2,0.6)`

Alternatively, a function `report_all` which uses three energy ranges (selected for `nicer`) by default:

```
help?> JAXTAM.report_all
  report_all(::Mission, ::DataFrameRow; e_ranges=[(0.2,12.0), (2.0,10.0), (0.2,2.0)], overwrite::Bool, nuke::Bool, update_masterpage::Bool)

  Calls report with three default energy ranges
```

### 6. Automated Reports

An `auto_report` function exists:

```
help?> JAXTAM.auto_report
  auto_report(::Mission; limit::Union{Bool,Int}, update::Bool, nuke::Bool)

  Calls auto_queue function to generate a queue of reports to make, the queue filters:

  * Public-only
  
  * Reportless
  
  * Not 'CAL' type observations
  
  * Error free

  Leaving only suitable observations to be analysed

  Calls 'report_all' using the queued observations

  Will continue to generate reports until the limit is reached (if there is one)
```

This will queue up some observations which meet the filter criterea mentioned above, then run `report_all` on them

# Appendix

Installing on Asimov should show:

```
srv01039:/home/robertr> cd Projects/
srv01039:/home/robertr/Projects> git clone https://github.com/RobertRosca/JAXTAM.jl
Cloning into 'JAXTAM.jl'...
remote: Enumerating objects: 169, done.
remote: Counting objects: 100% (169/169), done.
remote: Compressing objects: 100% (104/104), done.
remote: Total 2354 (delta 100), reused 125 (delta 63), pack-reused 2185
Receiving objects: 100% (2354/2354), 2.23 MiB | 2.22 MiB/s, done.
Resolving deltas: 100% (1512/1512), done.
Checking connectivity... done.
srv01039:/home/robertr/Projects> julia
[ Info: Running headless plots:  'ENV["GKSwstype"] = "100"'
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.1.0 (2019-01-21)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(v1.1) pkg> activate JAXTAM.jl/

(JAXTAM) pkg> instantiate
   Cloning default registries into `~/.julia`
   Cloning registry from "https://github.com/JuliaRegistries/General.git"
     Added registry `General` to `~/.julia/registries/General`
   Cloning git-repo `https://github.com/RobertRosca/Hyperscript.jl.git`
 Installed NaNMath ───────────────────── v0.3.2
 Installed PlotThemes ────────────────── v0.3.0
 Installed PenaltyFunctions ──────────── v0.1.2
 Installed SortingAlgorithms ─────────── v0.3.1
 Installed SpecialFunctions ──────────── v0.7.2
 Installed TranscodingStreams ────────── v0.9.3
 Installed LossFunctions ─────────────── v0.5.1
 Installed Colors ────────────────────── v0.9.5
 Installed LearnBase ─────────────────── v0.2.2
 Installed SweepOperator ─────────────── v0.2.0
 Installed VersionParsing ────────────── v1.1.3
 Installed BinaryProvider ────────────── v0.5.3
 Installed OnlineStats ───────────────── v0.20.3
 Installed JSON ──────────────────────── v0.20.0
 Installed Plots ─────────────────────── v0.23.0
 Installed FileIO ────────────────────── v1.0.5
 Installed QuadGK ────────────────────── v2.0.3
 Installed DataStreams ───────────────── v0.4.1
 Installed Compat ────────────────────── v2.1.0
 Installed Reexport ──────────────────── v0.2.0
 Installed PlotUtils ─────────────────── v0.5.5
 Installed DataStructures ────────────── v0.15.0
 Installed StaticArrays ──────────────── v0.10.3
 Installed LombScargle ───────────────── v0.4.0
 Installed BinDeps ───────────────────── v0.8.10
 Installed Conda ─────────────────────── v1.2.0
 Installed OnlineStatsBase ───────────── v0.9.3
 Installed TableTraits ───────────────── v0.4.1
 Installed Showoff ───────────────────── v0.2.1
 Installed RecipesBase ───────────────── v0.6.0
 Installed Parameters ────────────────── v0.10.3
 Installed NamedTuples ───────────────── v5.0.0
 Installed URIParser ─────────────────── v0.4.0
 Installed OrderedCollections ────────── v1.0.2
 Installed Arrow ─────────────────────── v0.2.3
 Installed FlatBuffers ───────────────── v0.5.3
 Installed ColorTypes ────────────────── v0.7.5
 Installed DataFrames ────────────────── v0.17.1
 Installed StatsBase ─────────────────── v0.27.0
 Installed Requires ──────────────────── v0.5.2
 Installed Measurements ──────────────── v2.0.0
 Installed Contour ───────────────────── v0.5.1
 Installed FFTW ──────────────────────── v0.2.4
 Installed CategoricalArrays ─────────── v0.5.2
 Installed CodecZlib ─────────────────── v0.5.2
 Installed FITSIO ────────────────────── v0.13.0
 Installed Feather ───────────────────── v0.5.1
 Installed GR ────────────────────────── v0.38.1
 Installed Tables ────────────────────── v0.1.18
 Installed AbstractFFTs ──────────────── v0.3.2
 Installed JLD2 ──────────────────────── v0.1.2
 Installed Polynomials ───────────────── v0.5.2
 Installed DSP ───────────────────────── v0.5.2
 Installed Missings ──────────────────── v0.4.0
 Installed FixedPointNumbers ─────────── v0.5.3
 Installed Measures ──────────────────── v0.3.0
 Installed WeakRefStrings ────────────── v0.5.8
 Installed IteratorInterfaceExtensions ─ v0.1.1
 Installed Calculus ──────────────────── v0.4.1
  Building SpecialFunctions → `~/.julia/packages/SpecialFunctions/fvheQ/deps/build.log`
  Building GR ──────────────→ `~/.julia/packages/GR/IVBgs/deps/build.log`
  Building Plots ───────────→ `~/.julia/packages/Plots/UQI78/deps/build.log`
  Building Conda ───────────→ `~/.julia/packages/Conda/CpuvI/deps/build.log`
  Building FFTW ────────────→ `~/.julia/packages/FFTW/p7sLQ/deps/build.log`
  Building CodecZlib ───────→ `~/.julia/packages/CodecZlib/9jDi1/deps/build.log`
  Building FITSIO ──────────→ `~/.julia/packages/FITSIO/2H5Bk/deps/build.log`

julia> using JAXTAM
[ Info: Precompiling JAXTAM [c0c225ea-a005-11e8-11c7-71dab99cc9f2]
[ Info: Precompiling DataFrames [a93c6f00-e57d-5684-b7b6-d8193f3e46c0]
[ Info: Precompiling FileIO [5789e2e9-d7fb-5bc7-8068-2c6fae9b9549]
[ Info: Precompiling JSON [682c06a0-de6a-54ab-a142-c8b1cf79cde6]
[ Info: Precompiling JLD2 [033835bb-8acc-5ee8-8aae-3f567f8a3819]
[ Info: Precompiling FITSIO [525bcba6-941b-5504-bd06-fd0dc1a4d2eb]
[ Info: Precompiling Arrow [69666777-d1a9-59fb-9406-91d4454c9d45]
[ Info: Precompiling Feather [becb17da-46f6-5d3c-ad1b-1c5fe96bc73c]
[ Info: Precompiling FFTW [7a1cc6ca-52ef-59f5-83cd-3a7055c09341]
[ Info: Precompiling OnlineStats [a15396b6-48d5-5d58-9928-6d29437db91e]
WARNING: Method definition std(OnlineStatsBase.OnlineStat{T} where T) in module OnlineStatsBase at /home/robertr/.julia/packages/OnlineStatsBase/x5KCe/src/OnlineStatsBase.jl:128 overwritten in module OnlineStats at /home/robertr/.julia/packages/OnlineStats/NseHX/src/utils.jl:34.
WARNING: Method definition #std(Any, typeof(Statistics.std), OnlineStatsBase.OnlineStat{T} where T) in module OnlineStatsBase overwritten in module OnlineStats.
[ Info: Precompiling LombScargle [fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9]
[ Info: Precompiling DSP [717857b8-e6f2-59f4-9121-6e50c889abd2]
[ Info: Precompiling Hyperscript [61f73626-e61e-46c8-aded-7ef67e964bac]
[ Info: Precompiling Measures [442fdcdd-2543-5da2-b0f3-8c86c306513e]
[ Info: Precompiling Plots [91a5bcdd-55d7-5caf-9e0b-520d859cae80
```