#!/usr/bin/env node
/**
 * update-repo-tree.js
 *
 * Generates a repository file-tree (similar to the Unix `tree` command) and
 * inserts it into README.md between the markers:
 *
 *   <!-- REPO-TREE-START -->
 *   <!-- REPO-TREE-END -->
 *
 * Usage (from repo root):
 *   node scripts/update-repo-tree.js
 */

'use strict';

const fs   = require('fs');
const path = require('path');

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
const REPO_ROOT  = path.resolve(__dirname, '..');
const README     = path.join(REPO_ROOT, 'README.md');
const MARKER_START = '<!-- REPO-TREE-START -->';
const MARKER_END   = '<!-- REPO-TREE-END -->';

// Directories/files to exclude entirely (by exact name).
const EXCLUDE_NAMES = new Set([
  '.git',
  'node_modules',
  'dist',
  'build',
  '.venv',
  '__pycache__',
  '.Rproj.user',
]);

// ---------------------------------------------------------------------------
// Tree generation
// ---------------------------------------------------------------------------

/**
 * Returns a sorted list of entries: directories first, then files,
 * each group sorted case-insensitively.
 */
function sortedEntries(dirPath) {
  let entries;
  try {
    entries = fs.readdirSync(dirPath, { withFileTypes: true });
  } catch (err) {
    console.warn(`WARNING: Could not read directory ${dirPath}: ${err.message}`);
    return [];
  }
  const dirs  = entries.filter(e => e.isDirectory() && !EXCLUDE_NAMES.has(e.name));
  const files = entries.filter(e => !e.isDirectory() && !EXCLUDE_NAMES.has(e.name));
  const cmp = (a, b) => a.name.toLowerCase().localeCompare(b.name.toLowerCase());
  dirs.sort(cmp);
  files.sort(cmp);
  return [...dirs, ...files];
}

/**
 * Recursively builds tree lines.
 *
 * @param {string}   dirPath  - absolute path of current directory
 * @param {string}   prefix   - indentation prefix for this level
 * @param {string[]} lines    - accumulator
 */
function buildTree(dirPath, prefix, lines) {
  const entries = sortedEntries(dirPath);
  entries.forEach((entry, idx) => {
    const isLast      = idx === entries.length - 1;
    const connector   = isLast ? '└── ' : '├── ';
    const childPrefix = isLast ? '    ' : '│   ';
    const suffix      = entry.isDirectory() ? '/' : '';
    lines.push(`${prefix}${connector}${entry.name}${suffix}`);
    if (entry.isDirectory()) {
      buildTree(path.join(dirPath, entry.name), prefix + childPrefix, lines);
    }
  });
}

function generateTree() {
  const lines = ['.'];
  buildTree(REPO_ROOT, '', lines);
  return lines.join('\n');
}

// ---------------------------------------------------------------------------
// README update
// ---------------------------------------------------------------------------

function updateReadme(tree) {
  const content = fs.readFileSync(README, 'utf8');
  const startIdx = content.indexOf(MARKER_START);
  const endIdx   = content.indexOf(MARKER_END);

  if (startIdx === -1 || endIdx === -1) {
    console.error(
      `ERROR: Could not find markers in README.md.\n` +
      `Expected:\n  ${MARKER_START}\n  ${MARKER_END}`
    );
    process.exit(1);
  }

  const before    = content.slice(0, startIdx + MARKER_START.length);
  const after     = content.slice(endIdx);
  const newContent = `${before}\n\`\`\`text\n${tree}\n\`\`\`\n${after}`;

  fs.writeFileSync(README, newContent, 'utf8');
  console.log('README.md updated successfully.');
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
const tree = generateTree();
updateReadme(tree);
