#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_8.pro
create_lutable_pressure_8

EOF
echo '------------  end  ------------'