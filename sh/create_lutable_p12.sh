#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_12.pro
create_lutable_pressure_12

EOF
echo '------------  end  ------------'