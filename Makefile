default: error

coverage:
	@coverage run --source pyschism -m nose --rednose --nologcapture --verbose && coverage report -m && coverage-badge -f -o tests/coverage.svg && rm -rf .coverage

badge:
	@$(shell git diff --exit-code -s tests/coverage.svg)
	@if ! [ $(.SHELLSTATUS) = 0 ]; then git add tests/coverage.svg && git commit -m "Updated coverage badge."; fi

error:
	@echo -e "ERROR: Argument required."
	@exit 2