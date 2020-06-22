#!/bin/bash

# fix these
declare errors="E301,E303,E304,E305,E306,E502,W291,W293,W391,E226,E225,E241,E231,"

# Be belligerent for these
declare belligerent="W292,E265,"

# Only warn for these
declare warnings="E101,E111,E112,E113,E121,E122,E125,E127,E128,E129,E131,E133,E201,E202,E203,E211,E222,E223,E224,E225,E227,E228,E242,E251,E261,E262,E271,E272,E402,E703,E711,E712,E713,E714,E721,W504,E302,YTT101,YTT102,YTT103,YTT201,YTT202,YTT203,YTT204,YTT301,YTT302,YTT303,STRFTIME001,STRFTIME002,PT001,PT002,PT003,PT004,PT005,PT006,PT007,PT008,PT009,PT010,PT011,PT012,PT013,PT014,PT015,PT016,PT017,PT018,PT019,PT020,PT021,RST201,RST202,RST203,RST204,RST205,RST206,RST207,RST208,RST210,RST211,RST212,RST213,RST214,RST215,RST216,RST217,RST218,RST219,RST299,RST301,RST302,RST303,RST304,RST305,RST306,RST399,RST401,RST499,RST900,RST901,RST902,RST903,Q000,Q001,Q002,Q003,"


if [ -z "$(git status --porcelain --untracked-files=no)" ] || [ "$1" == "-f" ]; then
  # Working directory clean

  echo "Running autopep8"

  autopep8 --in-place --select "$errors" -a --recursive pyms/
  autopep8 --in-place --select "$belligerent" -a -a -a -a -a --recursive pyms/

  autopep8 --in-place --select "$errors" -a --recursive tests/
  autopep8 --in-place --select "$belligerent" -a -a -a -a -a --recursive tests/

  echo "Running flake8"

    >&2 flake8 pyms/

    >&2 flake8 tests/

  exit 0

else
  # Uncommitted changes
  >&2 echo "git working directory is not clean"
  exit 1

fi
