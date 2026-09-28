"use strict";

const LIMITS = Object.freeze({ equationLength: 1000, formulaLength: 256, compounds: 32, elements: 64, nullity: 3, search: 12 });
const abs = value => value < 0n ? -value : value;
const gcd = (a, b) => { a = abs(a); b = abs(b); while (b !== 0n) [a, b] = [b, a % b]; return a || 1n; };
const lcm = (a, b) => abs(a * b) / gcd(a, b);

class Fraction {
  constructor(numerator = 0n, denominator = 1n) {
    if (denominator === 0n) throw new Error("Internal division by zero.");
    if (denominator < 0n) { numerator = -numerator; denominator = -denominator; }
    const factor = gcd(numerator, denominator);
    this.n = numerator / factor;
    this.d = denominator / factor;
  }
  add(other) { return new Fraction(this.n * other.d + other.n * this.d, this.d * other.d); }
  sub(other) { return new Fraction(this.n * other.d - other.n * this.d, this.d * other.d); }
  mul(other) { return new Fraction(this.n * other.n, this.d * other.d); }
  div(other) { return new Fraction(this.n * other.d, this.d * other.n); }
  neg() { return new Fraction(-this.n, this.d); }
  isZero() { return this.n === 0n; }
}

/** Parse a chemical formula into an element -> atom count Map. */
export function parseFormula(raw) {
  const formula = String(raw).trim();
  if (!formula) throw new Error("A chemical formula is missing.");
  if (formula.length > LIMITS.formulaLength) throw new Error("A chemical formula is too long.");
  if (!/^[A-Za-z0-9()[\]]+$/.test(formula)) throw new Error(`“${formula}” contains an unsupported character.`);
  let position = 0;
  const add = (counts, element, count) => counts.set(element, (counts.get(element) || 0n) + count);
  const readNumber = () => {
    const start = position;
    while (/\d/.test(formula[position] || "")) position++;
    if (start === position) return 1n;
    const value = BigInt(formula.slice(start, position));
    if (value < 1n) throw new Error(`Subscripts in “${formula}” must be positive whole numbers.`);
    return value;
  };
  const parseGroup = closing => {
    const counts = new Map();
    let found = false;
    while (position < formula.length && formula[position] !== closing) {
      const character = formula[position];
      if (character === "(" || character === "[") {
        const expected = character === "(" ? ")" : "]";
        position++;
        const nested = parseGroup(expected);
        if (formula[position] !== expected) throw new Error(`A closing ${expected} is missing in “${formula}”.`);
        position++;
        const multiplier = readNumber();
        nested.forEach((count, element) => add(counts, element, count * multiplier));
        found = true;
      } else if (character === ")" || character === "]") {
        throw new Error(`There is a mismatched bracket in “${formula}”.`);
      } else if (/[A-Z]/.test(character)) {
        let element = formula[position++];
        while (/[a-z]/.test(formula[position] || "")) element += formula[position++];
        add(counts, element, readNumber());
        found = true;
      } else if (/[a-z]/.test(character)) {
        throw new Error(`Element symbols must start with a capital letter in “${formula}”.`);
      } else {
        throw new Error(`A number in “${formula}” must follow an element or bracket.`);
      }
    }
    if (!found) throw new Error(closing ? `Empty brackets are not allowed in “${formula}”.` : `“${formula}” is not a valid formula.`);
    return counts;
  };
  const counts = parseGroup(null);
  if (position !== formula.length) throw new Error(`There is an unmatched bracket in “${formula}”.`);
  return counts;
}

/** Parse supported equation arrows and terms. Leading coefficients are accepted but do not constrain balancing. */
export function parseEquation(raw) {
  const text = String(raw).trim();
  if (!text) throw new Error("Enter a chemical equation to get started.");
  if (text.length > LIMITS.equationLength) throw new Error("The equation is too long.");
  const arrows = [...text.matchAll(/-->|->|→|⇌|=/g)];
  if (arrows.length === 0) throw new Error("Use ->, -->, =, →, or ⇌ to separate reactants and products.");
  if (arrows.length !== 1) throw new Error("Use exactly one arrow to separate reactants and products.");
  const arrow = arrows[0];
  const sides = [text.slice(0, arrow.index), text.slice(arrow.index + arrow[0].length)];
  if (!sides[0].trim() || !sides[1].trim()) throw new Error("Add at least one formula on each side of the arrow.");
  const parseSide = side => side.split("+").map(piece => {
    const term = piece.trim();
    if (!term) throw new Error("A formula is missing next to a plus sign.");
    const match = term.match(/^(\d+)?\s*([A-Za-z0-9()[\]]+)$/);
    if (!match) throw new Error(`“${term}” is not a valid chemical formula.`);
    if (match[1] && BigInt(match[1]) < 1n) throw new Error("Existing coefficients must be positive whole numbers.");
    return { formula: match[2], atoms: parseFormula(match[2]), inputCoefficient: match[1] ? BigInt(match[1]) : 1n };
  });
  const parsed = { reactants: parseSide(sides[0]), products: parseSide(sides[1]) };
  const count = parsed.reactants.length + parsed.products.length;
  if (count > LIMITS.compounds) throw new Error(`Equations may contain at most ${LIMITS.compounds} formulas.`);
  return parsed;
}

function nullSpace(matrix) {
  const reduced = matrix.map(row => row.map(value => new Fraction(value)));
  const rows = reduced.length;
  const columns = reduced[0].length;
  const pivots = [];
  let row = 0;
  for (let column = 0; column < columns && row < rows; column++) {
    let pivot = row;
    while (pivot < rows && reduced[pivot][column].isZero()) pivot++;
    if (pivot === rows) continue;
    [reduced[row], reduced[pivot]] = [reduced[pivot], reduced[row]];
    const divisor = reduced[row][column];
    reduced[row] = reduced[row].map(value => value.div(divisor));
    for (let other = 0; other < rows; other++) {
      if (other === row || reduced[other][column].isZero()) continue;
      const factor = reduced[other][column];
      reduced[other] = reduced[other].map((value, index) => value.sub(factor.mul(reduced[row][index])));
    }
    pivots.push(column);
    row++;
  }
  const free = Array.from({ length: columns }, (_, index) => index).filter(column => !pivots.includes(column));
  return free.map(freeColumn => {
    const vector = Array.from({ length: columns }, () => new Fraction());
    vector[freeColumn] = new Fraction(1n);
    pivots.forEach((pivotColumn, index) => { vector[pivotColumn] = reduced[index][freeColumn].neg(); });
    return vector;
  });
}

function toSmallestIntegers(vector) {
  let denominator = 1n;
  vector.forEach(value => { denominator = lcm(denominator, value.d); });
  let values = vector.map(value => value.n * (denominator / value.d));
  const factor = values.reduce((common, value) => gcd(common, value), 0n);
  values = values.map(value => value / factor);
  if (values.every(value => value < 0n)) values = values.map(value => -value);
  return values;
}

function compareSolutions(left, right) {
  const leftSum = left.reduce((sum, value) => sum + value, 0n);
  const rightSum = right.reduce((sum, value) => sum + value, 0n);
  if (leftSum !== rightSum) return leftSum < rightSum ? -1 : 1;
  for (let index = 0; index < left.length; index++) if (left[index] !== right[index]) return left[index] < right[index] ? -1 : 1;
  return 0;
}

function solveBasis(basis) {
  if (basis.length === 1) return toSmallestIntegers(basis[0]);
  if (basis.length > LIMITS.nullity) throw new Error("This equation has too many independent reactions to simplify automatically.");
  let best = null;
  const choices = Array(basis.length).fill(0n);
  const search = index => {
    if (index < choices.length) {
      for (let value = -LIMITS.search; value <= LIMITS.search; value++) {
        choices[index] = BigInt(value);
        search(index + 1);
      }
      return;
    }
    if (choices.every(value => value === 0n)) return;
    const vector = basis[0].map((_, compound) => basis.reduce(
      (sum, item, basisIndex) => sum.add(item[compound].mul(new Fraction(choices[basisIndex]))), new Fraction()
    ));
    const candidate = toSmallestIntegers(vector);
    if (candidate.every(value => value > 0n) && (!best || compareSolutions(candidate, best) < 0)) best = candidate;
  };
  search(0);
  if (!best) throw new Error("This equation cannot be balanced with positive coefficients.");
  return best;
}

export function isMassConserved(parsed, coefficients) {
  const leftCount = parsed.reactants.length;
  const elements = new Set([...parsed.reactants, ...parsed.products].flatMap(term => [...term.atoms.keys()]));
  return [...elements].every(element => {
    const left = parsed.reactants.reduce((sum, term, index) => sum + (term.atoms.get(element) || 0n) * coefficients[index], 0n);
    const right = parsed.products.reduce((sum, term, index) => sum + (term.atoms.get(element) || 0n) * coefficients[index + leftCount], 0n);
    return left === right;
  });
}

export function solveEquation(parsed) {
  const leftElements = new Set(parsed.reactants.flatMap(term => [...term.atoms.keys()]));
  const rightElements = new Set(parsed.products.flatMap(term => [...term.atoms.keys()]));
  if (leftElements.size !== rightElements.size || [...leftElements].some(element => !rightElements.has(element))) {
    throw new Error("Reactants and products must contain the same elements.");
  }
  if (leftElements.size > LIMITS.elements) throw new Error(`Equations may contain at most ${LIMITS.elements} elements.`);
  const elements = [...leftElements].sort();
  const compounds = [...parsed.reactants, ...parsed.products];
  const matrix = elements.map(element => compounds.map((term, index) =>
    (term.atoms.get(element) || 0n) * (index < parsed.reactants.length ? 1n : -1n)
  ));
  const basis = nullSpace(matrix);
  if (!basis.length) throw new Error("This equation cannot be balanced as written.");
  const coefficients = solveBasis(basis);
  if (coefficients.some(value => value <= 0n) || !isMassConserved(parsed, coefficients)) {
    throw new Error("This equation cannot be balanced with positive coefficients.");
  }
  return { coefficients, elements, matrix };
}

export function formatEquation(parsed, coefficients, arrow = "->") {
  const format = (term, index) => `${coefficients[index] === 1n ? "" : coefficients[index]}${term.formula}`;
  const left = parsed.reactants.map((term, index) => format(term, index)).join(" + ");
  const right = parsed.products.map((term, index) => format(term, index + parsed.reactants.length)).join(" + ");
  return `${left} ${arrow} ${right}`;
}

/** Pure entry point used by both the browser and unit tests. */
export function balanceEquation(equation) {
  const parsed = parseEquation(equation);
  const solution = solveEquation(parsed);
  return { ...solution, parsed, equation: formatEquation(parsed, solution.coefficients) };
}
