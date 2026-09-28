import assert from "node:assert/strict";
import { performance } from "node:perf_hooks";
import test from "node:test";

import { balanceEquation, isMassConserved, parseEquation, parseFormula } from "../src/lib/balancer.js";

const cases = [
  ["combustion", "CH4 + O2 -> CO2 + H2O", "CH4 + 2O2 -> CO2 + 2H2O"],
  ["synthesis", "N2 + H2 -> NH3", "N2 + 3H2 -> 2NH3"],
  ["decomposition", "KClO3 -> KCl + O2", "2KClO3 -> 2KCl + 3O2"],
  ["single replacement", "Fe + HCl -> FeCl2 + H2", "Fe + 2HCl -> FeCl2 + H2"],
  ["double replacement", "AgNO3 + NaCl -> AgCl + NaNO3", "AgNO3 + NaCl -> AgCl + NaNO3"],
  ["polyatomic group", "Al + CuSO4 -> Al2(SO4)3 + Cu", "2Al + 3CuSO4 -> Al2(SO4)3 + 3Cu"],
];

for (const [name, input, expected] of cases) {
  test(`balances ${name}`, () => {
    const result = balanceEquation(input);
    assert.equal(result.equation, expected);
    assert.equal(isMassConserved(result.parsed, result.coefficients), true);
  });
}

test("accepts every documented arrow", () => {
  for (const arrow of ["=", "->", "-->", "→", "⇌"]) {
    assert.equal(balanceEquation(`H2 + O2 ${arrow} H2O`).equation, "2H2 + O2 -> 2H2O");
  }
});

test("normalizes an already balanced equation and a trivial equation", () => {
  assert.equal(balanceEquation("C + O2 -> CO2").equation, "C + O2 -> CO2");
  assert.equal(balanceEquation("H2 -> H2").equation, "H2 -> H2");
  assert.equal(balanceEquation("2 H2 + 1 O2 -> 2 H2O").equation, "2H2 + O2 -> 2H2O");
});

test("parses nested parentheses and square brackets", () => {
  assert.deepEqual(Object.fromEntries(parseFormula("K4[ON(SO3)2]2")), { K: 4n, O: 14n, N: 2n, S: 4n });
});

test("builds a signed atom matrix", () => {
  const result = balanceEquation("H2 + O2 -> H2O");
  assert.deepEqual(result.elements, ["H", "O"]);
  assert.deepEqual(result.matrix, [[2n, 0n, -2n], [0n, 2n, -1n]]);
});

test("rejects mismatched element sets", () => {
  assert.throws(() => balanceEquation("H2 + O2 -> NaCl"), /same elements/);
});

test("rejects impossible positive balances", () => {
  assert.throws(() => balanceEquation("CO -> CO2 + CO3"), /positive coefficients/);
});

test("rejects incomplete and malformed equations", () => {
  const invalid = ["", "H2 + O2", "H2 + -> H2O", "H2 ->", "H2 -> -> H2O", "Mg(OH2 -> MgO", "2h -> H2", "0H2 -> H2", "H2 + O2 -> H2O +"];
  for (const equation of invalid) assert.throws(() => balanceEquation(equation), Error, equation);
});

test("guarded failures and the standard suite finish under 300ms", () => {
  const start = performance.now();
  for (let iteration = 0; iteration < 25; iteration++) {
    for (const [, equation] of cases) balanceEquation(equation);
    assert.throws(() => balanceEquation("H2 + O2 -> NaCl"));
  }
  assert.ok(performance.now() - start < 300, "balancing workload exceeded 300ms");
});

test("parser exposes reactants and products", () => {
  const parsed = parseEquation("CH4 + O2 --> CO2 + H2O");
  assert.deepEqual(parsed.reactants.map(term => term.formula), ["CH4", "O2"]);
  assert.deepEqual(parsed.products.map(term => term.formula), ["CO2", "H2O"]);
});
