#!/usr/bin/env bash
set -e

link=data_TF_japonica
expression=TF_japonica
water_class=all #CK
data_type=RNA
adj_type=adj #_drop
species=japonica
ID=japonica_test
#python main_graph.py name=GENELink_${water_class}_${expression}_${ID} dset=${dset} expression=${expression} model=GENELink batch_size=256 epochs=200 flag=True water_class=${water_class} data_type=${data_type} adj_type=${adj_type}
python main_graph.py name=test_${water_class}_${expression}_${ID} link=${link} expression=${expression} model=Graphormer2Link_gate_type batch_size=512 epochs=1 flag=True water_class=${water_class} data_type=${data_type} adj_type=${adj_type}





