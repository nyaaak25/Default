#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_7.pro
create_lutable_pressure_7

EOF
echo '------------  end  ------------'