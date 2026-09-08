# Mandatory Rules

## In all languages

- Prefer slightly verbose, self-explanatory code over terse code that needs
  comments to be understood.
- Keep comments to only what explains something non-obvious.
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

## Python

- In Python that uses ROOT, never use matplotlib. Look at nearby files for the
  established plotting approach, or ask which is preferred.

## Explanations

- For non-trivial changes, explain thoroughly what changed and why. Do not
  over-summarize or truncate the reasoning. Trivial edits can stay terse.
