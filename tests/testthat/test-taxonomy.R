test_that("ranks are assigned by prefix, not by position", {
  # Truncated at genus, the way Greengenes2 actually writes an unresolved
  # species. Position-based splitting would put the genus in `species`.
  tax <- ap_parse_taxonomy(
    "d__Bacteria; p__Bacillota_A_368345; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae; g__Blautia_A_141781"
  )
  expect_equal(tax$phylum, "Bacillota_A_368345")
  expect_equal(tax$genus, "Blautia_A_141781")
  expect_true(is.na(tax$species))
})

test_that("a full seven-rank string parses to seven ranks", {
  tax <- ap_parse_taxonomy(
    "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Phocaeicola_A; s__Phocaeicola_A vulgatus"
  )
  expect_equal(unname(unlist(tax)),
               c("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales",
                 "Bacteroidaceae", "Phocaeicola_A", "Phocaeicola_A vulgatus"))
})

test_that("a string missing intermediate ranks does not shift the deeper ones", {
  tax <- ap_parse_taxonomy("d__Bacteria; p__Bacteroidota; g__Bacteroides")
  expect_equal(tax$phylum, "Bacteroidota")
  expect_true(is.na(tax$class))
  expect_true(is.na(tax$order))
  expect_true(is.na(tax$family))
  expect_equal(tax$genus, "Bacteroides")
})

test_that("the Greengenes 13_8 k__ prefix is read as domain", {
  tax <- ap_parse_taxonomy("k__Bacteria; p__Firmicutes")
  expect_equal(tax$domain, "Bacteria")
  expect_equal(tax$phylum, "Firmicutes")
})

test_that("placeholders that name nothing become NA", {
  tax <- ap_parse_taxonomy(c(
    "d__Bacteria; p__Bacteroidota; g__uncultured",
    "d__Bacteria; p__Bacteroidota; g__",
    "d__Bacteria; p__Bacteroidota; g__metagenome",
    "Unassigned"
  ))
  expect_true(all(is.na(tax$genus)))
  expect_equal(tax$phylum, c("Bacteroidota", "Bacteroidota", "Bacteroidota", NA))
  expect_true(is.na(tax$domain[4]))
})

test_that("a genuine taxon whose name contains a placeholder word survives", {
  tax <- ap_parse_taxonomy("d__Bacteria; g__Unknownia")
  expect_equal(tax$genus, "Unknownia")
})

test_that("an empty string yields all NA rather than an error", {
  tax <- ap_parse_taxonomy(c("", "   "))
  expect_equal(nrow(tax), 2L)
  expect_true(all(is.na(unlist(tax))))
})

test_that("an unknown rank is refused by name", {
  expect_error(ap_parse_taxonomy("d__Bacteria", ranks = c("phylum", "kingdom")),
               "kingdom")
})

# --- labels ---

test_that("labels fall back up the hierarchy and say that they did", {
  tax <- ap_parse_taxonomy(c(
    "d__Bacteria; p__Bacteroidota; c__C; o__O; f__Bacteroidaceae; g__Bacteroides",
    "d__Bacteria; p__Bacteroidota; c__C; o__O; f__Lachnospiraceae",
    "Unassigned"
  ))
  lab <- ap_taxon_label(tax, "genus")
  expect_equal(lab[1], "g__Bacteroides")
  expect_equal(lab[2], "f__Lachnospiraceae (unassigned genus)")
  expect_equal(lab[3], "Unassigned")
})

test_that("fallback marking can be turned off but the name still comes from the right rank", {
  tax <- ap_parse_taxonomy("d__Bacteria; p__Bacteroidota; f__Lachnospiraceae")
  expect_equal(ap_taxon_label(tax, "genus", mark_fallback = FALSE),
               "f__Lachnospiraceae")
})

test_that("labelling at a rank absent from the table is refused", {
  tax <- ap_parse_taxonomy("d__Bacteria", ranks = c("domain", "phylum"))
  expect_error(ap_taxon_label(tax, "genus"), "not a column")
})
