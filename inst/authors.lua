-- Lua filter: render structured author/affiliation/orcid metadata
-- into LaTeX using the authblk package, suppressing Pandoc's default
-- \author{Name1 \and Name2} rendering.

function Meta(meta)
  if not meta["by-author"] then return meta end

  local lines = {}

  -- \author[affil numbers]{Name \orcidlink{id}}
  for _, author in ipairs(meta["by-author"]) do
    local name = pandoc.utils.stringify(author.name.literal)

    local aff_nums = {}
    if author.affiliations then
      for _, aff in ipairs(author.affiliations) do
        table.insert(aff_nums, tostring(aff.number))
      end
    end

    local orcid_str = ""
    if author.orcid then
      orcid_str = "\\,\\orcidlink{" .. pandoc.utils.stringify(author.orcid) .. "}"
    end

    table.insert(lines, string.format(
      "\\author[%s]{%s%s}",
      table.concat(aff_nums, ","), name, orcid_str
    ))
  end

  -- \affil[number]{Institution}
  if meta["affiliations"] then
    for _, aff in ipairs(meta["affiliations"]) do
      table.insert(lines, string.format(
        "\\affil[%s]{%s}",
        tostring(aff.number),
        pandoc.utils.stringify(aff.name)
      ))
    end
  end

  -- neutralize any subsequent \author{} calls from the Pandoc template
  table.insert(lines, "\\renewcommand*{\\author}[2][]{}")

  -- inject into header-includes
  local raw = pandoc.RawBlock("latex", table.concat(lines, "\n"))
  if not meta["header-includes"] then
    meta["header-includes"] = pandoc.MetaList({ pandoc.MetaBlocks({ raw }) })
  else
    table.insert(meta["header-includes"], pandoc.MetaBlocks({ raw }))
  end

  -- suppress Pandoc's default \author{Name1 \and Name2} in the preamble
  -- (Quarto normalizes 'author' into 'authors', so both must be cleared)
  meta.author = nil
  meta.authors = nil

  return meta
end
