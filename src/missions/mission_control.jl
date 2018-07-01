"""
    MissionDefinition

Type holding critical mission variables: mission name,
HEASARC mastertable url, and mission folder path
"""
mutable struct MissionDefinition
    name::String
    url::String
    path::String
    path_obs::Function
    path_cl::Function
    path_uf::Function
end