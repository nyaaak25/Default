#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_13.pro
create_lutable_pressure_13

EOF
echo '------------  end  ------------'