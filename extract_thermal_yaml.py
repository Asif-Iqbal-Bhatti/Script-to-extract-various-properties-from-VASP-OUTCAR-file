#!/usr/bin/env python3

import yaml

with open("thermal_properties.yaml") as f:
	data = yaml.load(f, Loader=yaml.FullLoader)


# isinstance checks whether the value is a list and if list
# we count the number of elements it contains.
# dict = {key, value} pair

count = sum(len(data[x]) for x in data if isinstance(data[x], list))
print(f"# of elements detected in the list: {count}")


print ("{:1s} {:14.12s} {:14.12s} {}".format(' ','Temperature', 'Entropy', 'Free_energy') )
for i in range(count):
	temperature=data['thermal_properties'][i]['temperature']
	free_energy=data['thermal_properties'][i]['free_energy']
	entropy=data['thermal_properties'][i]['entropy']
	heat_capacity=data['thermal_properties'][i]['heat_capacity']
	energy=data['thermal_properties'][i]['energy']
	print ("{:10.3f} {:15.8f}  {:15.8f}".format(temperature, entropy, free_energy) )