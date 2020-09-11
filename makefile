##############################################################################
#
#   Makefile for the documentation of program PB3D (Peeling Ballooning in 3D)
#   Toon Weyens
#
##############################################################################

##############################################################################
#   Variables
##############################################################################
export PB3D_dir = /home/tweyens/Simulations/PB3D
PB3D_version := $(shell grep 'prog_version =' $(PB3D_dir)/Modules/num_vars.f90 | cut --complement -d = -f 1 | sed -e 's/^ *//g' | cut -d'_' -f1)
MIN_NM_X := $(shell grep 'min_nm_X =' $(PB3D_dir)/Modules/X_vars.f90 | cut --complement -d = -f 1 | sed -e 's/^ *//g' | cut -d' ' -f1)
DOXY_version := $(shell doxygen --version)

##############################################################################
#   Rules
##############################################################################
clean:
	rm -f -r Doxygen/*

all: doc doc_latex

git-push:
	git push --follow-tags

info:
	@echo 'PROJECT_NUMBER = $(PB3D_version)'
	@echo 'ALIASES += min_nm_X="$(MIN_NM_X)"'
	@echo 'ALIASES += doxy_version="$(DOXY_version)"'

doc: 
	@set -e
	@rm -f temp_user_vars
	@echo 'PROJECT_NUMBER = [$(PB3D_version)]' >> temp_user_vars
	@echo 'ALIASES += min_nm_X="$(MIN_NM_X)"' >> temp_user_vars
	@echo 'ALIASES += doxy_version="$(DOXY_version)"' >> temp_user_vars
	( cat Doxyfile temp_user_vars ) | doxygen -
	@rm -f temp_user_vars
	@./Scripts/clean_html_and_latex.sh

doc_latex: doc
	cd ./Doxygen/latex && lualatex refman.tex && bibtex refman && lualatex refman.tex && lualatex refman.tex && mv refman.pdf ../../PB3D_manual.pdf
	@echo "\n Created file 'PB3D_manual.pdf'."
	@echo "\n Warnings in refman.log file in './Doxygen/latex':\n"
	@echo "START OF WARNINGS"
	@cd ./Doxygen/latex && grep -i 'warning' refman.log
	@echo "END OF WARNINGS"
