#!/usr/bin/env Rscript
# render_docs.R
# Converts Markdown developer guides into standalone, beautifully styled HTML pages
# for publication via GitHub Pages (docs/ directory).

docs_dir <- "docs"
if (!dir.exists(docs_dir)) dir.create(docs_dir, recursive = TRUE)

files_to_render <- c(
  "DEVELOPER.md"    = "Developer Guide — MafA Discovery",
  "shinyapp.md"     = "Legacy Standalone Prototypes — MafA Discovery",
  "publishapp.md"   = "Publishing & Deployment Guide — MafA Discovery",
  "redesign.md"     = "UI Redesign & Reactivity — MafA Discovery",
  "qtlanalysis.md"  = "F2 Glycemic QTL Integration — MafA Discovery"
)

# Convert Markdown to HTML fragment
md_to_html <- function(md_content) {
  if (requireNamespace("commonmark", quietly = TRUE)) {
    commonmark::markdown_html(
      md_content,
      extensions = c("table", "strikethrough", "autolink", "tagfilter")
    )
  } else if (requireNamespace("markdown", quietly = TRUE)) {
    markdown::markdownToHTML(text = md_content, fragment.only = TRUE)
  } else {
    # Basic fallback
    paste0("<pre>", htmltools::htmlEscape(md_content), "</pre>")
  }
}

html_template <- function(title, body_html, current_file) {
  nav_links <- c(
    "DEVELOPER.html"   = "Developer Guide",
    "shinyapp.html"    = "Legacy Prototypes",
    "publishapp.html"  = "Publishing & WebAssembly",
    "redesign.html"    = "UI Redesign",
    "qtlanalysis.html" = "QTL Analysis"
  )
  
  nav_items <- vapply(names(nav_links), function(f) {
    active_cls <- if (f == current_file) ' class="active"' else ''
    sprintf('<li><a href="%s"%s>%s</a></li>', f, active_cls, nav_links[[f]])
  }, character(1))
  nav_html <- paste(nav_items, collapse = "\n          ")
  
  sprintf('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>%s</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
  <style>
    :root {
      --primary: #1976d2;
      --primary-hover: #115293;
      --text: #24292f;
      --text-muted: #57606a;
      --bg: #ffffff;
      --bg-alt: #f6f8fa;
      --border: #d0d7de;
      --code-bg: #f6f8fa;
      --accent: #6f42c1;
    }
    * { box-sizing: border-box; }
    body {
      font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      line-height: 1.65;
      color: var(--text);
      background-color: #fafbfc;
      margin: 0;
      padding: 0;
    }
    header.site-header {
      background: #ffffff;
      border-bottom: 1px solid var(--border);
      position: sticky;
      top: 0;
      z-index: 100;
      box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    }
    .header-inner {
      max-width: 1060px;
      margin: 0 auto;
      padding: 12px 24px;
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 12px;
    }
    .header-brand {
      font-size: 1.15em;
      font-weight: 700;
      color: var(--primary);
      text-decoration: none;
      display: flex;
      align-items: center;
      gap: 8px;
    }
    .header-brand:hover { color: var(--primary-hover); }
    nav.doc-nav ul {
      list-style: none;
      margin: 0;
      padding: 0;
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    nav.doc-nav a {
      color: var(--text-muted);
      text-decoration: none;
      font-size: 0.9em;
      font-weight: 500;
      padding: 6px 12px;
      border-radius: 6px;
      transition: all 0.15s ease;
    }
    nav.doc-nav a:hover {
      background-color: var(--bg-alt);
      color: var(--primary);
    }
    nav.doc-nav a.active {
      background-color: #e3f2fd;
      color: var(--primary);
      font-weight: 600;
    }
    .app-link {
      background-color: #2da44e !important;
      color: #ffffff !important;
      font-weight: 600 !important;
    }
    .app-link:hover {
      background-color: #2c974b !important;
    }
    main.doc-container {
      max-width: 1060px;
      margin: 32px auto 60px auto;
      padding: 0 24px;
    }
    article.markdown-body {
      background: #ffffff;
      padding: 48px;
      border: 1px solid var(--border);
      border-radius: 8px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.03);
    }
    @media (max-width: 768px) {
      article.markdown-body { padding: 24px 16px; }
      .header-inner { padding: 12px 16px; }
    }
    h1, h2, h3, h4, h5, h6 {
      color: #1f2328;
      font-weight: 600;
      margin-top: 28px;
      margin-bottom: 14px;
      line-height: 1.3;
    }
    h1 {
      font-size: 2.1em;
      border-bottom: 1px solid var(--border);
      padding-bottom: 12px;
      margin-top: 0;
    }
    h2 {
      font-size: 1.5em;
      border-bottom: 1px solid var(--border);
      padding-bottom: 8px;
    }
    h3 { font-size: 1.25em; }
    p, ul, ol { margin-top: 0; margin-bottom: 16px; }
    ul, ol { padding-left: 28px; }
    li { margin-bottom: 6px; }
    a { color: var(--primary); text-decoration: underline; text-underline-offset: 2px; }
    a:hover { color: var(--primary-hover); }
    code {
      font-family: "JetBrains Mono", monospace;
      font-size: 0.88em;
      background-color: var(--code-bg);
      padding: 3px 6px;
      border-radius: 4px;
      border: 1px solid #e1e4e8;
    }
    pre {
      background-color: var(--code-bg);
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 16px;
      overflow-x: auto;
      font-size: 0.9em;
      line-height: 1.5;
    }
    pre code {
      background: none;
      padding: 0;
      border: none;
    }
    table {
      width: 100%%;
      border-collapse: collapse;
      margin: 20px 0;
      font-size: 0.94em;
    }
    th, td {
      border: 1px solid var(--border);
      padding: 10px 14px;
      text-align: left;
    }
    th {
      background-color: var(--bg-alt);
      font-weight: 600;
    }
    tr:nth-child(even) td { background-color: #fcfcfd; }
    blockquote {
      border-left: 4px solid var(--primary);
      margin: 16px 0;
      padding: 8px 16px;
      background-color: #f0f7ff;
      color: #333;
    }
    hr {
      border: 0;
      border-top: 1px solid var(--border);
      margin: 32px 0;
    }
    footer.site-footer {
      text-align: center;
      padding: 24px;
      color: var(--text-muted);
      font-size: 0.88em;
      border-top: 1px solid var(--border);
      background: #ffffff;
    }
  </style>
</head>
<body>
  <header class="site-header">
    <div class="header-inner">
      <a href="../index.html" class="header-brand">🧬 MafA Discovery</a>
      <nav class="doc-nav">
        <ul>
          <li><a href="../index.html" class="app-link">🚀 Open App</a></li>
          %s
          <li><a href="https://github.com/byandell/MafADiscovery" target="_blank" rel="noopener">GitHub ↗</a></li>
        </ul>
      </nav>
    </div>
  </header>
  
  <main class="doc-container">
    <article class="markdown-body">
%s
    </article>
  </main>

  <footer class="site-footer">
    <p>MafA Discovery Integrated Genomic Explorer &bull; Vanderbilt University / University of Wisconsin-Madison</p>
  </footer>
</body>
</html>', title, nav_html, body_html)
}

cat("Rendering markdown documentation to HTML...\n")

for (md_file in names(files_to_render)) {
  if (!file.exists(md_file)) {
    warning("File not found: ", md_file)
    next
  }
  
  title <- files_to_render[[md_file]]
  target_html <- file.path(docs_dir, sub("\\.md$", ".html", md_file))
  current_file <- basename(target_html)
  
  cat(" -", md_file, "->", target_html, "\n")
  md_content <- paste(readLines(md_file, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  
  # Adjust internal relative links for HTML output:
  # Replace .md links with .html links
  md_content <- gsub("\\]\\(DEVELOPER\\.md\\)", "](DEVELOPER.html)", md_content)
  md_content <- gsub("\\]\\(shinyapp\\.md\\)", "](shinyapp.html)", md_content)
  md_content <- gsub("\\]\\(publishapp\\.md\\)", "](publishapp.html)", md_content)
  md_content <- gsub("\\]\\(redesign\\.md\\)", "](redesign.html)", md_content)
  md_content <- gsub("\\]\\(qtlanalysis\\.md\\)", "](qtlanalysis.html)", md_content)
  
  body_html <- md_to_html(md_content)
  full_html <- html_template(title, body_html, current_file)
  
  writeLines(full_html, target_html, useBytes = TRUE)
}

# Ensure .nojekyll exists
writeLines("", file.path(docs_dir, ".nojekyll"))

cat("Documentation successfully rendered to docs/*.html\n")
