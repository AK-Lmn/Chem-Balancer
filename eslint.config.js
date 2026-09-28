import js from "@eslint/js";

export default [
  {
    ignores: ["node_modules/**"],
  },
  js.configs.recommended,
  {
    files: ["script.js"],
    languageOptions: {
      globals: {
        document: "readonly",
        localStorage: "readonly",
        navigator: "readonly",
        setTimeout: "readonly",
      },
    },
  },
  {
    files: ["test/**/*.js"],
    languageOptions: {
      globals: {
        BigInt: "readonly",
      },
    },
  },
];
