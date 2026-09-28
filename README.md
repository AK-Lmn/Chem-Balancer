# Chem-Balancer

Chem-Balancer is a browser-based chemical equation balancer. The interface calls a standalone, pure JavaScript engine in `src/lib/balancer.js`; it has no DOM or browser dependencies and can be tested directly in Node.js.

## Balancing algorithm

1. The parser separates both sides of the reaction and recursively counts atoms in each formula, including nested parentheses and square brackets.
2. The engine creates a signed element-by-compound matrix: reactant counts are positive and product counts are negative.
3. Exact rational Gaussian elimination finds the matrix null space, avoiding floating-point rounding.
4. Rational values are converted to their smallest positive integer ratio and independently checked for mass conservation before being returned.

Malformed formulas, mismatched element sets, equations with no positive solution, and excessively large or underdetermined inputs return a clear error. Fixed bounds on input size and multi-solution search prevent runaway work.

## Supported input

- Reaction separators: `=`, `->`, `-->`, `→`, and `⇌`
- Compounds separated with `+`
- Standard case-sensitive element symbols and positive integer subscripts
- Nested groups such as `Al2(SO4)3` and `K4[ON(SO3)2]2`
- Optional existing whole-number coefficients, which are normalized to the smallest ratio

Example: `CH4 + O2 -> CO2 + H2O` becomes `CH4 + 2O2 -> CO2 + 2H2O`.

Ionic charges, electron notation, phases, hydrates, and isotope notation are not currently supported.

## Development and tests

Node.js 20.19 or newer is required. The balancing engine has no runtime dependencies; ESLint is used for development checks.

```sh
npm ci
npm run check
npm test
```

The deterministic `node:test` suite covers combustion, synthesis, decomposition, replacement reactions, polyatomic groups, all supported arrows, trivial/already-balanced reactions, atom matrices, syntax failures, impossible equations, mass conservation, and a sub-300ms workload guard.
