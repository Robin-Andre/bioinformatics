#!/bin/bash

find ./practical-2021/test/res/reference_results -iname 'RAxML_RF-Distances.0' -execdir mv -i '{}' distances \;
find ./practical-2021/test/res/reference_results -iname 'RAxML_info.0' -execdir mv -i '{}' info \;
