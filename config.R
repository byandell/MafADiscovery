# config.R
# MafA Discovery — Central Repository Configuration
# Single Source of Truth for repository ownership and URLs.
#
# When transferring this repository to a new owner or organization:
# 1. Update `owner` below (e.g., owner = "new-owner").
# 2. Run `Rscript render_docs.R` to sync all Markdown links and recompile docs/*.html.
#
# Architectural guide & checklist: https://byandell.github.io/DocumentTools/prompts/move_repo.html

GITHUB_CONFIG <- list(
  owner  = "byandell",          # GitHub user or organization
  repo   = "MafADiscovery",     # Repository name
  branch = "main"               # Default production branch
)

# Derived URLs
GITHUB_CONFIG$pages_url <- sprintf("https://%s.github.io/%s", GITHUB_CONFIG$owner, GITHUB_CONFIG$repo)
GITHUB_CONFIG$repo_url  <- sprintf("https://github.com/%s/%s", GITHUB_CONFIG$owner, GITHUB_CONFIG$repo)
GITHUB_CONFIG$docs_url  <- sprintf("%s/docs", GITHUB_CONFIG$pages_url)
GITHUB_CONFIG$blob_url  <- sprintf("%s/blob/%s", GITHUB_CONFIG$repo_url, GITHUB_CONFIG$branch)

#' Synchronize absolute GitHub and GitHub Pages links across repository Markdown files
#' Scans documentation files and ensures all absolute GitHub and Pages links reflect the current owner.
sync_markdown_links <- function(config = GITHUB_CONFIG, dry_run = FALSE) {
  md_files <- c(
    "README.md",
    "DEVELOPER.md",
    list.files("guides", pattern = "\\.md$", full.names = TRUE)
  )
  md_files <- md_files[file.exists(md_files)]
  
  updated_count <- 0
  
  for (f in md_files) {
    content <- readLines(f, encoding = "UTF-8", warn = FALSE)
    orig_text <- paste(content, collapse = "\n")
    new_text <- orig_text
    
    # Replace GitHub Pages URLs: https://<any-owner>.github.io/<repo>/...
    pages_regex <- sprintf("https://[a-zA-Z0-9_-]+\\.github\\.io/%s", config$repo)
    new_text <- gsub(pages_regex, config$pages_url, new_text)
    
    # Replace GitHub repository URLs: https://github.com/<any-owner>/<repo>
    repo_regex <- sprintf("https://github\\.com/[a-zA-Z0-9_-]+/%s", config$repo)
    new_text <- gsub(repo_regex, config$repo_url, new_text)
    
    # Replace markdown repo slugs: [<any-owner>/<repo>]
    slug_regex <- sprintf("\\[[a-zA-Z0-9_-]+/%s\\]", config$repo)
    new_text <- gsub(slug_regex, sprintf("[%s/%s]", config$owner, config$repo), new_text)
    
    if (new_text != orig_text) {
      if (!dry_run) {
        writeLines(strsplit(new_text, "\n")[[1]], f, useBytes = TRUE)
      }
      updated_count <- updated_count + 1
      cat(sprintf("  [%s] %s\n", if (dry_run) "Would update" else "Updated", f))
    }
  }
  
  invisible(updated_count)
}
