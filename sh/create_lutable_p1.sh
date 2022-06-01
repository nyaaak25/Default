#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_1.pro
create_lutable_pressure_1

EOF
echo '------------  end  ------------'