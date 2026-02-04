#!/bin/bash
#
# # Check if a filename was provided
# if [ -z "$1" ]; then
#     echo "Usage: $0 <filename.R>"
#     exit 1
# fi
#
# FILE_PATH=$1
# PROJ=${pwd}
FILE_PATH="${CAP_PROJECT_PATH}/src/12_drugtarget_analysis.R"

# Run lintr and styler inside the container
singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_models_0.1.1.sif R -e \
                 "styler::style_file('$FILE_PATH'); lintr::lint('$FILE_PATH')"

