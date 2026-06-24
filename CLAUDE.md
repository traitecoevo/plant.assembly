# plant.assembly — developer guide for Claude

R package for **assembling ecological communities with the `plant` model**: it
repeatedly invades empty/occupied trait space, solves each community to
demographic equilibrium, computes invasion-fitness landscapes, and finds
evolutionary attractors (selection gradients, singular strategies).

It is a thin evolutionary-assembly layer **on top of `plant`**. `plant` provides
the demographic model (the SCM — Size- and Patch-structured Cohort Model);
`plant.assembly` orchestrates many `plant` runs.

## Status (June 2026)

The package was mid-way through a refactor that **ports equilibrium/demography
code out of `plant` and into this package**, while `plant` itself went through
breaking interface changes. As of now:

- `devtools::load_all()` **succeeds**.
- The **core path works and is verified**: `community_start` → `community_add` →
  `community_demography()` (default `equilibrium_iteration` solver) →
  `community_selection_gradient()`. At equilibrium resident fitness ≈ 0 and the
  selection gradient is finite — matches `vignettes/solving_attractors.Rmd`.
- **Viable-bounds path works** and is reimplemented on the community machinery:
  `community_viable_fitness_1D(community)` / `community_viable_bounds()` use an
  empty community's `fitness_function` as the fundamental-fitness function (no
  dependency on plant's removed `fundamental_fitness()`/`viable_fitness()`).
- See **Known issues** below for what is *not* yet working.

### Important: plant removed the whole fitness/equilibrium subsystem (#388)

plant NEWS (#388): "All fitness/equilibrium functionality was removed from plant;
it now lives in `plant.assembly`." Removed from plant: `fitness_landscape()`,
`solve_max_fitness()`, `viable_fitness()`, `fundamental_fitness()`,
`equilibrium_birth_rate()`, `positive_1d()`, `max_growth_rate()`, etc. The port
is *exactly* the job of bringing these into this package. The intended approach
is to **reimplement on the community machinery** (run the SCM via the community's
`fitness_function` / `community_demography`) rather than copy plant's old code
verbatim. Done so far: equilibrium (→ `community_demography`), viable bounds
(→ `community_viable_fitness_1D`), `positive_1d`/`positive_1d_bracket` (pure
numeric helpers, in `R/community_fitness_viable.R`), and the fitness-max
functions `max_fitness()` / `max_growth_rate()` (`R/community_fitness_solve_max.R`,
operating on `community$fitness_function`). The package no longer references any
removed/undefined plant symbol (verified with `codetools::findGlobals`).

`HANDOFF.md` is an ephemeral note from the previous session (not committed);
this file supersedes it for durable guidance.

## Dev commands

```r
devtools::load_all(".")     # load package (works)
devtools::document()        # regenerate man/ + NAMESPACE from roxygen
devtools::test()            # run testthat suite (see baseline below)
```

There is a manual smoke script mirroring the core vignette path; the canonical
end-to-end exercises live in `vignettes/solving_attractors.Rmd`,
`vignettes/assembly_fitmax.Rmd`, `vignettes/assembly_stochastic.Rmd`.

## Architecture

### The `community` object (`R/community.R`)

A list with class `community`. Built by `community_start(bounds, ...)`, then
mutated through a pipeline. Key fields:

- `bounds`, `trait_names`, `traits` (matrix), `birth_rate` (per resident)
- `demography_control` — how to solve demographic equilibrium (see controls)
- `model_support` — the bridge to `plant`: `list(p = <plant Parameters>,
  plant_control = <plant control()>)`. Also caches `node_schedule_times` /
  `node_schedule_ode_times` between solves.
- `fitness_control`, `fitness_function`, `fitness_points`, `resident_fitness`,
  `selection_gradient` — populated by the fitness/gradient functions.

The typical pipeline (see `solving_attractors.Rmd`):

```r
community_start(bounds, model_support = list(p = ..., plant_control = ...)) |>
  community_add(trait_matrix(x, "lma"), birth_rate = 200) |>
  community_demography() |>            # solve to equilibrium birth rates
  community_selection_gradient()       # fitness gradient at residents
```

### plant bridge (`R/community_plant.R`)

All direct `plant` coupling lives here. The package-internal generic names
(`community_parameters`, `community_make_demography_runner`,
`community_demography_runner_cleanup`, `community_viable_bounds`,
`community_check_for_inviable_strategies`) are **aliased** at the bottom of the
file to their `plant_community_*` implementations. This indirection is so the
plant-specific layer could in principle be swapped.

- `plant_default_assembly_pars(hmat, max_patch_lifetime, fixed_RA)` — base FF16
  `Parameters`. **Important:** it regenerates `node_schedule_times_default` to
  match `max_patch_lifetime` (otherwise the default schedule spans the original
  ~105 and the SCM errors with "time_max must be greater than current time").
- `plant_default_assembly_control(...)` — wrapper around `plant::control()`.
- `plant_community_parameters(community)` — builds a `plant` `Parameters` from
  the community (strategies via `plant::strategy_list`, restores cached schedule
  times).
- `plant_community_make_demography_runner(community)` — returns a closure
  `runner(birth_rates) -> offspring_production` that runs the SCM once. It
  stashes `last_schedule_times` / `last_offspring_production` / `history` in its
  environment for the cleanup step. Resets the integration schedule to defaults
  when birth rates jump by more than
  `demography_control$equilibrium_large_birth_rate_change`.
- `plant_community_demography_runner_cleanup(community, runner, converged)` —
  reads the runner's environment and writes the equilibrium `birth_rate`,
  schedule times, and `converged`/`progress` attrs back onto the community.

### Demography / equilibrium (`R/community_demography.R`)

`community_demography(community)` dispatches on
`demography_control$equilibrium_solver_name`:

- `single_step` — one SCM step, no iteration.
- `equilibrium_iteration` — **default, working**. Fixed-point iteration of
  incoming→outgoing offspring until `equilibrium_eps` is reached.
- `equilibrium_solve_nleqslv` / `equilibrium_solve_dfsane` — root-finding via
  `util_nlsolve` (`R/util_nlsolve.R`). **Partially ported / unverified.**
- `equilibrium_hybrid` — iterate then root-find. **Broken** (see Known issues).

After solving, `plant_community_update_fitness_function()` builds the mutant
invasion-fitness closure on the community.

### Fitness, gradients, assembly

- `R/community_fitness_landscape.R`, `community_fitness_viable.R`,
  `community_fitness_solve_max.R` — invasion-fitness landscapes (some use
  `mlr3`/Gaussian-process surrogates) and viable trait bounds.
- `R/solve_attractors.R` — `community_selection_gradient()`,
  `community_solve_singularity_1D()`. Finite-difference gradients use the
  internal `gradient_points()`/`gradient_extrapolate()` in `R/util_gradient.R`.
- `R/assembler.R` — `assembler_start`/`assembler_run`/`assembler_control` drive
  full assembly (births → demography → deaths) over many steps.
- `R/births*.R`, `R/deaths.R` — add/remove strategies (maximum-fitness or
  stochastic births; inviable-strategy removal).

### Control objects — keep them straight

| Object | Built by | Purpose |
|---|---|---|
| `demography_control` | `demographic_step_control()` | Equilibrium solving: `equilibrium_solver_name`, `equilibrium_eps`, `equilibrium_nsteps`, `equilibrium_large_birth_rate_change`, `equilibrium_extinct_birth_rate`, etc. Lives at `community$demography_control`. |
| `plant_control` | `plant_default_assembly_control()` / `plant::control()` | Passed straight to `run_scm()`. A plant `Control` S4 object — **cannot** hold extra fields, so all equilibrium params go in `demography_control`. Lives at `community$model_support$plant_control`. |
| `fitness_control` | (list) | How fitness landscapes are sampled (`method`, `n_evals`, …). |
| `assembler_control` | `assembler_control()` | Assembly loop: birth/death type, tolerances. |

Note: plant renamed `seed_rain` → `birth_rate`/`offspring`; control field names
here follow current plant terminology.

## plant dependency and interface

- Built against the post-#459 `plant` (installed `2.0.0.9001`). Baseline recorded
  in `.plant-interface-version`.
- Source checkouts of plant live alongside this repo, e.g.
  `../plant-dev1` (see its `NEWS.md` "Breaking changes" for old→new maps) and
  `../plant-master`.
- Key migrated calls: `scm_base_control()`→`control()`;
  `run_scm_collect(p)`→`run_scm(p, collect=TRUE)`; `build_schedule(p, ctrl)` +
  `attr(., "offspring_production")` → `run_scm(p, ctrl, refine_schedule=TRUE)`
  then `scm$parameters` / `scm$offspring_production`;
  `p$node_schedule_ode_times` (Parameters field) → `p$ode_times`.
- `run_scm(..., collect=FALSE)` returns the **SCM object** (`scm$parameters`,
  `scm$offspring_production`, `scm$net_reproduction_ratios`, `scm$run_mutant`);
  with `collect=TRUE` it returns a tidied list.
- After future plant updates, re-run the `plant-update-interface` skill (in
  `../plant-dev1/.claude/skills/`) — it reads plant's NEWS and updates
  `.plant-interface-version`.

## Test baseline

`devtools::test()` is **green: 79 pass, 0 fail, 0 skip, 0 warn**.

- `test-community.R` (new) — covers the community interface: `trait_matrix`,
  `bounds`, `demographic_step_control`, `community_start/add/drop`,
  `length.community`, the `max_patch_lifetime` schedule regression, and
  integration tests for `community_demography` (empty + single resident,
  reference birth rate ≈ 112.66) and `community_selection_gradient`.
- `test-support-fitness.R` — `positive_1d`, `bounds`/`check_bounds`/`check_point`,
  and `community_viable_fitness_1D` (viable interval ≈ [0.0131, 1.644] for lma).

The five old test files that called plant's removed API were resolved
(drop-or-implement, not skipped): `positive_1d` was reimplemented so its test was
kept; the fundamental/viable tests were rewritten against the community machinery;
`test-x_equilibriumR.R` (superseded by `test-community.R`), `test-support-assembly.R`
(`assembly_parameters`), `test-trait_fitness.R` (dead empty loop), and the
duplicate `test-fitness-support.R` were deleted.

## Known issues / TODO

- **2D maximum-fitness births (`find_max_fitness_2d`) are unverified** — the path
  is now self-contained (uses `sys$fitness_function`) but relies on
  `fitness_slopes`/`maximize_logspace` and is commented as "not very well tested".
  Only 1D fitmax births are exercised.
- **`equilibrium_hybrid` solver is broken**: it treats the iteration result as a
  plant `Parameters` (`eq_solution$strategies`) but it is now a `community`, and
  it passes a `ctrl=` arg that `demography_solve_equilibrium_solve(community,
  solver)` no longer accepts. Needs reworking for the community object.
- **`equilibrium_solve_*` (nleqslv/dfsane) are unverified** end-to-end on the
  community object — only `equilibrium_iteration` has been exercised.
- **`community_plot_fitness_landscape`** is not working.
- `fitness_control` has no default (it is `NULL` out of `community_start()`); the
  fitness-landscape functions that read it should supply sensible defaults.
- `plant_community_check_for_inviable_strategies` still has a TODO to drop its
  direct plant dependency and reuse the community fitness functions.
- Two malformed `@param` roxygen tags in `R/community_plots.R` (lines 64, 103)
  warn on `document()`.
- Nothing in the working tree is committed yet — ~70 dirty files mixing the port
  WIP and the plant-interface migration. Stage deliberately.

## Issue & project-board conventions

Development across `plant`, `plant.assembly`, and `overstorey` is tracked on a
shared [project board](https://github.com/orgs/traitecoevo/projects/5). New issues
are auto-added to the board with status **Backlog** by a workflow, so you do not
need to set status manually.

When opening an issue (including whenever the user asks you to create one), always:

- **Set exactly one type label.** Only three labels exist in these repos — do not
  invent new ones:
  - `bug` — an existing feature not functioning as intended
  - `task` — a discrete task needed for a feature (the default for normal work)
  - `epic` — a new feature or capability, usually an umbrella over several tasks
- **Prefix the title with a theme tag** in square brackets so the board sorts
  cleanly. Assembly work is usually `[evol assembly]`; reuse another existing
  theme where it fits, or fall back to `[other]`:

  | Tag | Scope |
  |---|---|
  | `[evol assembly]` | Evolutionary assembly linking plant to plant.assembly |
  | `[TF24 hydraulics]` | Hydraulics component of the TF24 strategy |
  | `[TF24 allometry]` | Flexible allometry for the TF24 model |
  | `[TF24 nsc]` | Non-structural carbohydrate storage in TF24 |
  | `[acclimation]` | Acclimation of leaf and other traits |
  | `[simplify interface]` | Consistent interface to the plant & plant.assembly models |
  | `[Env drivers]` | Driving the model with environmental drivers |
  | `[speed]` | Performance — making the model run faster |
  | `[patch variations]` | Multiple patch setups (multi-patch, stochastic metapopulation, continuous patch) |
  | `[forecasting]` | Enabling forecasting with the plant model |
  | `[documentation]` | Documenting model capabilities (any of the three repos) |
  | `[other]` | Anything not covered above |

  A title may carry more than one tag when it genuinely spans themes.

Create issues with `gh issue create -R traitecoevo/plant.assembly
--title "[evol assembly] …" --label task` (swap in `bug`/`epic` as appropriate).
