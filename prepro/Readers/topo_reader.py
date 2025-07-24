def lookvar(topo_file):
    # etopo5.nc from Jeroen's easygrid
    if "etopo5" in topo_file.split("/")[-1].lower():
        dico = {"lon": "topo_lon", "lat": "topo_lat", "topo": "topo", "zaxis": "up"}
    # etopo2.nc from Romstools
    elif "etopo2" in topo_file.split("/")[-1].lower():
        dico = {"lon": "lon", "lat": "lat", "topo": "topo", "zaxis": "up"}
    elif "etopo1" in topo_file.split("/")[-1].lower():
        dico = {"lon": "x", "lat": "y", "topo": "z", "zaxis": "up"}
    # srtm file
    elif "srtm30" in topo_file.split("/")[-2].lower():
        dico = {
            "lon": "longitude",
            "lat": "latitude",
            "topo": "elevation",
            "zaxis": "up",
            "srtm": True,
        }
    # homonym file :
    # for French atlantic coast data dowloaded from
    #  https://diffusion.shom.fr/donnees/bathymerie/mnt-facade-atl-homonim.html
    elif "homonim" in topo_file.split("/")[-1].lower():
        dico = {"lon": "x", "lat": "y", "topo": "z", "zaxis": "up"}
    # gebco file
    elif "gebco" in topo_file.split("/")[-1].lower():
        dico = {"lon": "lon", "lat": "lat", "topo": "topo", "zaxis": "up"}

    return dico
