#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_14.pro
create_lutable_pressure_14

EOF
echo '------------  end  ------------'