import pypact as pp

filename = "fispii.out"

with pp.Reader(filename) as output:
    for t in output.inventory_data:
        print(t.irradiation_time)
        print(t.currenttime)
        print(t.cooling_time)
