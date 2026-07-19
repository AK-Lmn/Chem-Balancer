"use strict";

const input = document.querySelector("#equation");
const errorBox = document.querySelector("#error-message");
const results = document.querySelector("#results");
const balancedEquation = document.querySelector("#balanced-equation");
const coefficientsBox = document.querySelector("#coefficients");
const atomTable = document.querySelector("#atom-table");
const HISTORY_KEY = "chemkit-balancer-history";
let lastPlainResult = "";

// A tiny exact fraction type. BigInt keeps Gaussian elimination free of
// floating-point rounding, which is essential for reliable coefficients.
const abs = n => n < 0n ? -n : n;
function gcd(a, b) { a = abs(a); b = abs(b); while (b) [a, b] = [b, a % b]; return a || 1n; }
function lcm(a, b) { return abs(a * b) / gcd(a, b); }
class Fraction {
  constructor(n = 0n, d = 1n) {
    if (d === 0n) throw new Error("Internal division by zero.");
    if (d < 0n) { n = -n; d = -d; }
    const common = gcd(n, d); this.n = n / common; this.d = d / common;
  }
  add(o) { return new Fraction(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { return new Fraction(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { return new Fraction(this.n * o.n, this.d * o.d); }
  div(o) { return new Fraction(this.n * o.d, this.d * o.n); }
  neg() { return new Fraction(-this.n, this.d); }
  isZero() { return this.n === 0n; }
}

// Recursive formula parser. parseGroup stops at ')' and its caller applies
// the following multiplier, allowing nested groups such as K4(ON(SO3)2)2.
function parseFormula(raw) {
  const formula = raw.trim();
  if (!formula) throw new Error("A chemical formula is missing.");
  if (!/^[A-Za-z0-9()]+$/.test(formula)) throw new Error(`“${formula}” contains an unsupported character.`);
  let pos = 0;
  const readNumber = () => {
    const start = pos;
    while (/\d/.test(formula[pos] || "")) pos++;
    if (start === pos) return 1n;
    const value = BigInt(formula.slice(start, pos));
    if (value < 1n) throw new Error(`Subscripts in “${formula}” must be positive whole numbers.`);
    return value;
  };
  const add = (map, element, count) => map.set(element, (map.get(element) || 0n) + count);
  const parseGroup = (insideParentheses = false) => {
    const counts = new Map();
    let found = false;
    while (pos < formula.length && formula[pos] !== ")") {
      if (formula[pos] === "(") {
        pos++; const nested = parseGroup(true);
        if (formula[pos] !== ")") throw new Error(`A closing parenthesis is missing in “${formula}”.`);
        pos++; const multiplier = readNumber();
        nested.forEach((count, element) => add(counts, element, count * multiplier)); found = true;
      } else if (/[A-Z]/.test(formula[pos])) {
        let element = formula[pos++];
        while (/[a-z]/.test(formula[pos] || "")) element += formula[pos++];
        const multiplier = readNumber(); add(counts, element, multiplier); found = true;
      } else if (formula[pos] === ")") break;
      else if (/[a-z]/.test(formula[pos])) throw new Error(`Element symbols must start with a capital letter in “${formula}”.`);
      else if (/\d/.test(formula[pos])) throw new Error(`A number in “${formula}” must follow an element or parenthesis.`);
      else throw new Error(`Could not read “${formula}” near character ${pos + 1}.`);
    }
    if (!found) throw new Error(insideParentheses ? `Empty parentheses are not allowed in “${formula}”.` : `“${formula}” is not a valid formula.`);
    return counts;
  };
  const counts = parseGroup();
  if (pos < formula.length) throw new Error(`There is an unmatched parenthesis in “${formula}”.`);
  return counts;
}

function parseEquation(text) {
  const normalized = text.trim().replace(/→/g, "->");
  const arrows = normalized.match(/->|=/g);
  if (!arrows) throw new Error("Use -> or = to separate reactants and products.");
  if (arrows.length !== 1) throw new Error("Use exactly one arrow to separate reactants and products.");
  const sides = normalized.split(/->|=/);
  if (!sides[0].trim() || !sides[1].trim()) throw new Error("Add at least one formula on each side of the arrow.");
  const parseSide = side => side.split("+").map(piece => {
    const term = piece.trim();
    if (!term) throw new Error("A formula is missing next to a plus sign.");
    const match = term.match(/^(\d+)?\s*([A-Za-z0-9()]+)$/);
    if (!match) throw new Error(`“${term}” is not a valid chemical formula.`);
    const formula = match[2];
    return { formula, atoms: parseFormula(formula) }; // Input coefficients do not constrain the new balance.
  });
  return { reactants: parseSide(sides[0]), products: parseSide(sides[1]) };
}

// Reduced row-echelon form yields an exact null-space basis for A·x = 0.
function nullSpace(matrix) {
  const a = matrix.map(row => row.map(n => new Fraction(n)));
  const rows = a.length, cols = a[0].length, pivots = [];
  let r = 0;
  for (let c = 0; c < cols && r < rows; c++) {
    let pivot = r; while (pivot < rows && a[pivot][c].isZero()) pivot++;
    if (pivot === rows) continue;
    [a[r], a[pivot]] = [a[pivot], a[r]];
    const divisor = a[r][c]; a[r] = a[r].map(v => v.div(divisor));
    for (let i = 0; i < rows; i++) if (i !== r && !a[i][c].isZero()) {
      const factor = a[i][c]; a[i] = a[i].map((v, j) => v.sub(factor.mul(a[r][j])));
    }
    pivots.push(c); r++;
  }
  const free = [...Array(cols).keys()].filter(c => !pivots.includes(c));
  return free.map(freeCol => {
    const vector = Array.from({ length: cols }, () => new Fraction());
    vector[freeCol] = new Fraction(1n);
    pivots.forEach((pivotCol, row) => { vector[pivotCol] = a[row][freeCol].neg(); });
    return vector;
  });
}

function toSmallestIntegers(vector) {
  let denominator = 1n;
  vector.forEach(v => { denominator = lcm(denominator, v.d); });
  let values = vector.map(v => v.n * (denominator / v.d));
  let common = values.reduce((g, n) => gcd(g, n), 0n); values = values.map(n => n / common);
  if (values.every(n => n < 0n)) values = values.map(n => -n);
  return values;
}

function balance(parsed) {
  const leftElements = new Set(parsed.reactants.flatMap(t => [...t.atoms.keys()]));
  const rightElements = new Set(parsed.products.flatMap(t => [...t.atoms.keys()]));
  if (leftElements.size !== rightElements.size || [...leftElements].some(e => !rightElements.has(e)))
    throw new Error("Reactants and products must contain the same elements.");
  const elements = [...leftElements].sort();
  const compounds = [...parsed.reactants, ...parsed.products];
  const matrix = elements.map(element => compounds.map((term, i) => (term.atoms.get(element) || 0n) * (i < parsed.reactants.length ? 1n : -1n)));
  const basis = nullSpace(matrix);
  if (!basis.length) throw new Error("This equation cannot be balanced as written.");
  if (basis.length === 1) {
    const answer = toSmallestIntegers(basis[0]);
    if (answer.some(n => n <= 0n)) throw new Error("This equation cannot be balanced with positive coefficients.");
    return { coefficients: answer, elements };
  }
  // Underdetermined equations have several valid balances. Search small positive
  // free-variable combinations and choose the solution with the smallest sum.
  let best = null;
  const tryValues = (index, choices) => {
    if (index < basis.length) { for (let n = 1n; n <= 12n; n++) tryValues(index + 1, [...choices, n]); return; }
    const vector = basis[0].map((_, i) => basis.reduce((sum, b, j) => sum.add(b[i].mul(new Fraction(choices[j]))), new Fraction()));
    const ints = toSmallestIntegers(vector);
    if (ints.every(n => n > 0n) && (!best || ints.reduce((a,b)=>a+b,0n) < best.reduce((a,b)=>a+b,0n))) best = ints;
  };
  if (basis.length <= 4) tryValues(0, []);
  if (!best) throw new Error("This equation has multiple independent reactions and could not be simplified automatically.");
  return { coefficients: best, elements };
}

function subscriptFormula(formula) { return formula.replace(/\d+/g, n => [...n].map(d => "₀₁₂₃₄₅₆₇₈₉"[Number(d)]).join("")); }
function formatEquation(parsed, coeffs, html = false) {
  const terms = [...parsed.reactants, ...parsed.products];
  const format = (term, i) => `${coeffs[i] === 1n ? "" : coeffs[i].toString()}${html ? subscriptFormula(term.formula) : term.formula}`;
  const left = parsed.reactants.map((t, i) => format(t, i)).join(" + ");
  const right = parsed.products.map((t, i) => format(t, i + parsed.reactants.length)).join(" + ");
  return `${left} ${html ? "→" : "->"} ${right}`;
}

function render(parsed, solution) {
  const terms = [...parsed.reactants, ...parsed.products];
  lastPlainResult = formatEquation(parsed, solution.coefficients);
  balancedEquation.textContent = formatEquation(parsed, solution.coefficients, true);
  coefficientsBox.innerHTML = solution.coefficients.map((c, i) => `<div class="coefficient">${c}<small>${subscriptFormula(terms[i].formula)}</small></div>`).join("");
  atomTable.innerHTML = solution.elements.map(element => {
    const total = (side, offset) => side.reduce((sum, term, i) => sum + (term.atoms.get(element) || 0n) * solution.coefficients[i + offset], 0n);
    const left = total(parsed.reactants, 0), right = total(parsed.products, parsed.reactants.length);
    return `<tr><td>${element}</td><td>${left}</td><td>${right}</td><td>✓ Equal</td></tr>`;
  }).join("");
  results.hidden = false; results.scrollIntoView({ behavior: "smooth", block: "center" });
  saveHistory(lastPlainResult); renderHistory();
}

function showError(message) { errorBox.textContent = message; errorBox.hidden = false; results.hidden = true; }
function runBalancer() {
  errorBox.hidden = true;
  if (!input.value.trim()) return showError("Enter a chemical equation to get started.");
  try { const parsed = parseEquation(input.value); render(parsed, balance(parsed)); }
  catch (error) { showError(error.message || "We could not balance that equation. Please check the formulas and try again."); }
}

function getHistory() { try { return JSON.parse(localStorage.getItem(HISTORY_KEY)) || []; } catch { return []; } }
function saveHistory(equation) {
  const items = getHistory().filter(item => item.equation !== equation);
  items.unshift({ equation, time: Date.now() }); localStorage.setItem(HISTORY_KEY, JSON.stringify(items.slice(0, 6)));
}
function renderHistory() {
  const items = getHistory(), list = document.querySelector("#history-list"), empty = document.querySelector("#empty-history");
  empty.hidden = items.length > 0;
  list.innerHTML = items.map((item, i) => `<button class="history-item" data-history="${i}"><code>${subscriptFormula(item.equation).replace("->", "→")}</code><time>${new Date(item.time).toLocaleDateString()}</time></button>`).join("");
  list.querySelectorAll("[data-history]").forEach(button => button.addEventListener("click", () => { input.value = items[Number(button.dataset.history)].equation; runBalancer(); }));
}

document.querySelector("#balance-button").addEventListener("click", runBalancer);
input.addEventListener("keydown", event => { if (event.key === "Enter") runBalancer(); });
document.querySelector("#clear-button").addEventListener("click", () => { input.value = ""; errorBox.hidden = true; results.hidden = true; input.focus(); });
document.querySelectorAll(".example").forEach(button => button.addEventListener("click", () => { input.value = button.dataset.equation; runBalancer(); }));
document.querySelector("#copy-button").addEventListener("click", async event => { await navigator.clipboard.writeText(lastPlainResult); event.currentTarget.textContent = "Copied!"; setTimeout(() => event.currentTarget.textContent = "Copy equation", 1300); });
document.querySelector("#clear-history").addEventListener("click", () => { localStorage.removeItem(HISTORY_KEY); renderHistory(); });
renderHistory();
