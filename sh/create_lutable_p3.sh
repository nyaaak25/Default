#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_3.pro
create_lutable_pressure_3

EOF
echo '------------  end  ------------'