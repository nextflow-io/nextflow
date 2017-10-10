
#
# AWS Batch tests
# 
set -e 
$NXF_CMD run hello -c awsbatch.config -profile hello
$NXF_CMD run rnatoy -c awsbatch.config -profile rnatoy
$NXF_CMD run https://github.com/CRG-CNAG/CalliNGS-NF -c awsbatch.config -profile callings
