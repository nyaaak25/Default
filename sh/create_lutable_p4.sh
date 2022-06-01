#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_4.pro
create_lutable_pressure_4

EOF
echo '------------  end  ------------'