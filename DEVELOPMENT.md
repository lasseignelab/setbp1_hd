Development Setup
===============
This document provides recommendations for setting up a development environment
to facilitate reproducible research.
# Github Rulesets
Github rulesets define rules for interacting with a repository.  They are
intended to improve code quality through better code review and merging
practices.

A Github ruleset can be setup by importing a ruleset JSON file.
To setup a ruleset follow these steps:
* Navigate in your browser to the Github repository for the project created
with the `cap new` command.
* Navigate to the file `.github/settings/default_ruleset.sh` in Github and
download it using the "Download raw file" icon.
* Navigate to the Github repository settings tab.
* Expand the Rules section.
* Select Rulesets.
* Click the New Ruleset button.
* Select Import a ruleset.
* Choose the `default_ruleset.json` file downloaded previously.
* Scroll to the bottom of the screen and click the Create button to save
the ruleset.
