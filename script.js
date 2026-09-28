"use strict";

import { balanceEquation, formatEquation } from "./src/lib/balancer.js";

const input = document.querySelector("#equation");
const errorBox = document.querySelector("#error-message");
const results = document.querySelector("#results");
const balancedEquation = document.querySelector("#balanced-equation");
const coefficientsBox = document.querySelector("#coefficients");
const atomTable = document.querySelector("#atom-table");
const HISTORY_KEY = "chemkit-balancer-history";
let lastPlainResult = "";

function subscriptFormula(formula) { return formula.replace(/\d+/g, n => [...n].map(d => "₀₁₂₃₄₅₆₇₈₉"[Number(d)]).join("")); }

function render(parsed, solution) {
  const terms = [...parsed.reactants, ...parsed.products];
  lastPlainResult = formatEquation(parsed, solution.coefficients);
  balancedEquation.textContent = subscriptFormula(formatEquation(parsed, solution.coefficients, "→"));
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
  try { const solution = balanceEquation(input.value); render(solution.parsed, solution); }
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
document.querySelector("#copy-button").addEventListener("click", async event => {
  const button = event.currentTarget;
  try {
    await navigator.clipboard.writeText(lastPlainResult);
    button.textContent = "Copied!";
  } catch {
    button.textContent = "Copy failed";
  }
  setTimeout(() => { button.textContent = "Copy equation"; }, 1300);
});
document.querySelector("#clear-history").addEventListener("click", () => { localStorage.removeItem(HISTORY_KEY); renderHistory(); });
renderHistory();
