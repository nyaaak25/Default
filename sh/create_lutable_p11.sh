#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_11.pro
create_lutable_pressure_11

EOF
echo '------------  end  ------------'