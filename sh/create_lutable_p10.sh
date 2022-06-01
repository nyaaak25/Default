#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_10.pro
create_lutable_pressure_10

EOF
echo '------------  end  ------------'