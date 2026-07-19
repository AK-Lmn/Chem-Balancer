# Chemical Equation Balancer — ChemKit

A standalone, beginner-friendly chemical equation balancer with ChemKit's dark, glassy science theme. It parses formulas, calculates exact smallest whole-number coefficients, and displays a per-element atom count check. Recent results are stored locally in the browser—there is no backend or database.

## Run locally

No installation or build step is required. Open `index.html` in a modern browser, or serve the folder with any static server, for example:

```sh
python -m http.server 8000
```

Then open `http://localhost:8000`.

## Supported syntax

- Separate reactants and products with `->`, `=`, or `→`.
- Separate compounds with `+`; surrounding spaces are optional.
- Use standard, case-sensitive element symbols such as `H`, `O`, `Fe`, and `Cl`.
- Use whole-number subscripts: `H2O`, `CO2`, `C6H12O6`.
- Use parentheses and nested groups: `Ca(OH)2`, `Al2(SO4)3`.
- Leading coefficients are accepted but recalculated to the smallest valid ratio.

## Examples

- `H2 + O2 -> H2O`
- `Fe + O2 -> Fe2O3`
- `C3H8 + O2 -> CO2 + H2O`
- `Al + HCl -> AlCl3 + H2`
- `Ca(OH)2 + HCl -> CaCl2 + H2O`

## Known limitations

- Ionic charges, electrons, hydrate-dot notation, phases such as `(aq)`, and isotope notation are not supported.
- Element symbols are parsed syntactically; the app does not validate them against the periodic table.
- Highly underdetermined equations with more than four independent solution variables require the equation to be split into individual reactions.
- History is browser-specific and is removed if site storage is cleared.
