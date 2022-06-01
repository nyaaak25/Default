#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_2.pro
create_lutable_pressure_2

EOF
echo '------------  end  ------------'