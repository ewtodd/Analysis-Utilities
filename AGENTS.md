# Mandatory Rules

## In all languages

- Prefer slightly verbose, self-explanatory code over terse code that needs
  comments to be understood.
- Keep comments to only what explains something non-obvious. This governs
  inline commentary inside function bodies. Doxygen blocks on public
  declarations are exempt and have their own rules under "API documentation
  (Doxygen)" below.
- Never embed a literal `\n` inside a string or print argument. A line break
  is always its own explicit statement. In C++/ROOT, use
  `std::cout << ... << std::endl;`. In Python, split output into separate
  `print()` calls, and use a bare `print()` for a blank line rather than
  appending `\n`.

## Nix

- Always use flakes and flake-based commands (`nix run`, `nix shell`, etc.).
  Never use the old `nix-shell` approach.
- If you are confused, stop and ask for help. This is especially critical in
  Nix.
- Follow the existing style of the surrounding modules.

## C++ / ROOT

- Use ROOT data types, and pick the *correct* one for the actual need rather
  than defaulting blindly: `Int_t` for ordinary ints, `Long64_t` for entry
  counts and large/64-bit values, `Double_t` for floating point, `TString`
  for string convenience, and so on. Match the width and signedness the code
  actually requires.
- Do not use modern C++ features: no `auto`, no smart pointers, no
  range-based (`for (x : c)`) iteration. Use explicit types and classic
  indexed/iterator loops.
- Lambdas are permitted where they are short, local, and improve readability
  over the alternatives: sort comparators defined next to the std::sort call,
  thread workers whose explicit parameter lists would be longer and harder to
  scan than an inline capture. A lambda that spans more than about 5 lines or
  captures by reference outside an immediately obvious scope (e.g. stored in a
  std::function returned from the function) should still be a named function.
  When in doubt, write a named function.
- In performance-critical code, always gate logging behind a compile- or
  run-time toggle so it can be disabled. The `std::endl` flush is therefore
  never a concern on hot paths.

## API documentation (Doxygen)

- Doxygen blocks on public declarations are API reference material, not inline
  commentary. They are read in the generated HTML with no source next to them,
  so the "only what is non-obvious" rule does not apply: document the obvious
  parameters too. A reader who can already tell from the source what `sigma`
  means is not the reader these are written for.
- Public declarations in `include/*.hpp` — classes, structs, enums, free
  functions, and public member functions — carry a Doxygen block. Coverage is
  complete and `WARN_IF_UNDOCUMENTED` is on, so adding an undocumented public
  declaration fails `nix build .#docs`. Private members are not extracted and
  do not need blocks.
- Document a function at one declaration only. A forward declaration and its
  real declaration are the same entity to doxygen, so two blocks get merged and
  the `@param` list is reported as doubled — see the note above
  `LaunchInteractiveFitEditor` in `FittingUtils.hpp`.
- Document at minimum what the entity does, every `@param` (with units and
  valid ranges where they matter), `@return`, and `@throws` where it applies.
  Add `@note` or `@warning` for ownership, thread-safety, and ROOT global-state
  side effects — `gDirectory`/`gPad` changes, and which `TFile` owns a returned
  histogram.
- Use `///` for one-liners and `/** ... */` for longer blocks.
  `JAVADOC_AUTOBRIEF` is on, so the first sentence becomes the brief and
  everything after it the detailed description.
- Document the declaration in the header, never the definition in `src/`.
  Comments inside `src/*.cpp` remain ordinary implementation notes and stay
  subject to the "only what is non-obvious" rule.
- Build with `make docs` or `cmake --build build --target docs`. Output lands
  in `docs/` (gitignored), warnings in `docs/doxygen-warnings.log`. That log is
  currently empty and must stay that way — `nix build .#docs` fails outright if
  it is not, and that derivation is what the live site is served from.

## Documentation layout

- `doc/pages/*.md` is the narrative documentation — guides, worked examples,
  build and deployment notes. It is what the README used to hold. Prose that
  explains *how to use* the library goes here, not in the README.
- `doc/pages/index.md` is the site's main page and the only place the guide list
  lives; a new page must be added there with `@subpage <id>` or it will not
  appear in the navigation.
- Page ids are set with `{#id}` after the title. Do not reuse a directory name
  from `INPUT` as a page id — doxygen resolves the directory first and the
  `@subpage` silently points at a source listing instead of the page.
- `*.md` must stay in `FILE_PATTERNS`. Without it the directory scan skips every
  page, and the only symptom is unresolved `@subpage` warnings.
- `README.md` is a short pointer to the published site for people browsing
  GitHub. It is deliberately not in doxygen's `INPUT`; do not grow it back into
  a second copy of the documentation.
- The site's theme (the Kanagawa palette) lives in the website repository
  (github:ewtodd/website), not here. Change colours there.

## Python

- In Python that uses ROOT, never use matplotlib. Look at nearby files for the
  established plotting approach, or ask which is preferred.

## Explanations

- For non-trivial changes, explain thoroughly what changed and why. Do not
  over-summarize or truncate the reasoning. Trivial edits can stay terse.

## Publishing

- docs.ethanwtodd.com/analysis-utilities is the website repository's `docs`
  package (~/Software/website), which takes this repository as a flake input;
  the Caddy vhost on `nu` serves it. There is no CI publish step.
- To publish a documentation change: push here, run `nix flake update
  analysis-utilities` in ~/Software/website and push that, then `nix flake
  update website` in /etc/nixos and deploy nu (`colmena apply --on nu`).
- Doxygen emits only relative links, which is what allows the site to be served
  under a path prefix. Do not introduce root-absolute URLs into `doc/pages`, or
  the prefixed deployment breaks while local previews still work.

## Review checklists

- Every change to core analysis logic is human-reviewed before commit
  (README). Core: fitting models (`FittingUtils`, `RooFitUtils`,
  `RooFitPhotopeakPdfs` and the PDFs' `doEval`/integrals, fit math in the
  `Interactive*FitEditor` classes), signal processing
  (`WaveformProcessingUtils`), result extraction (`FitResult`,
  `Save*/Load*InteractiveParams`). `BinaryUtils`, `InitUtils`, `IOUtils`,
  `PlottingUtils`, `gpu/*.cu` dispatch, build files and docs are tooling. A
  small fix to fit logic is core; a `gpu/*.cu`-only perf change is tooling
  unless it changes the math.
- For a core change, write `checks/REVIEW-<slug>-<YYYY-MM-DD>.md` (tracked),
  following the newest file there. No code: `path:line-range`, what the lines
  do, and a `Check:` asking ewt to confirm something specific, not a
  restatement. Line numbers from the working tree (`grep -n`), not a diff;
  ascending within a file; core first, then tooling, then "checked and
  excluded", nits last. Leave boxes `- [ ]`.
- Agent commits are the ones authored by son-of-anton-bot (the history carries
  no `Co-Authored-By` trailers).

## Commits

- A Conventional Commits summary line with scope (`feat(PlottingUtils): …`,
  `fix(RooFitUtils): …`, `docs: …`), a blank line, then a body that explains
  what changed and why.
- Agent commits: `--author="son-of-anton-bot <307402699+son-of-anton-bot@users.noreply.github.com>"`.
