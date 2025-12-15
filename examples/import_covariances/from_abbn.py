import sauna

# From ABBN
covariances = sauna.Covariances()
covariances.library = 'ABBN'
covariances.group_structure = sauna.GROUP_STRUCTURES['ICSBEP']
covariances.from_abbns('./Abbn_Mf81') # Provide the path to the folder with covariances

print(covariances)
print(covariances.covariances[0])