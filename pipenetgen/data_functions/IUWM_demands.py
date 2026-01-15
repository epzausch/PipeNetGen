import osmnx as ox
import pandas as pd
import geopandas as gpd

import geopandas as gpd
import censusdis.data as ced
from censusdis.datasets import ACS5
from censusdis.states import NAMES_FROM_IDS


def read_IUWM(filename):
    df = pd.read_csv(filename)
    return df

def fix_IUWM(df):
    df['indoor_outdoor'] = df['indoor.demand'] + df['outdoor.demand']
    df['day'] = 1
    df['year_month'] = pd.to_datetime(df[['year','month','day']])
    return df

def find_max_month(df):
    gb = df.groupby("year_month")['indoor_outdoor'].max()
    month = gb.idxmax()
    value = gb.max()
    demands = df[df['year_month']==month].copy()
    demands['days'] = demands['year_month'].dt.days_in_month 
    return demands

def find_demand(df,peaking_factor=3):
    GPM_to_LPS = 0.0630902
    df.loc[:,'LPS'] = df['indoor_outdoor']/df['days']/24/60 * GPM_to_LPS
    df.loc[:,"peak_LPS"] = df['LPS'] * peaking_factor
    return df


def get_census_info(CENSUS_API_KEY,YEAR=2014,DATASET=ACS5,STATE="California",CONCEPT="Total Population",GROUP="B01003"):
    IDS_FROM_NAMES = {val:key for key,val in NAMES_FROM_IDS.items()}
    STATE = IDS_FROM_NAMES[STATE]

    df_acs5 = ced.download(
        DATASET,
        YEAR,
        leaves_of_group=GROUP,
        state=STATE,
        block_group="*",
        api_key=CENSUS_API_KEY,
    )

    return df_acs5

def match_geoids(demand_df,census_df):
    # Split the geoids so we can match to IUWM input

    split_geoids = pd.DataFrame(census_df.GEO_ID.str.split('US').tolist(),
                                    columns = ['num','geoid'])
    # print(split_geoids)
    # Add to df
    census_df['geoid'] = split_geoids["geoid"]

    # Make types match between IUWM and this df
    demand_df.loc[:, 'geoid'] = demand_df['geoid'].astype(str)

    match_df = census_df.loc[census_df['geoid'].isin(demand_df['geoid'])]
    if len(match_df) == 0:
        print("Something went wrong, trying US0 for split")
        split_geoids = pd.DataFrame(census_df.GEO_ID.str.split('US0').tolist(),
                                    columns = ['num','geoid'])
        # print(split_geoids)
        # Add to df
        census_df['geoid'] = split_geoids["geoid"]

        # Make types match between IUWM and this df
        demand_df['geoid'] = demand_df['geoid'].astype(str)

        match_df = census_df.loc[census_df['geoid'].isin(demand_df['geoid'])]
    
    data_df = match_df.merge(demand_df, on="geoid") # Merge data

    return data_df

def add_geometry_to_census(df,YEAR=2014,STATE=None):
    gdf_di = ced.add_inferred_geography(df, YEAR) # Get geoms
    return gdf_di


def match_node_to_census(G,df):
    for node in G.nodes:
        G.nodes[node]['count'] = 1 # For counting how many nodes per index
        G.nodes[node]['name'] = node # Common identifier that won't change
    
    nodes,edges = ox.graph_to_gdfs(G)
    nodes = nodes.to_crs(crs=df.crs)
    # print(df)
    # Number of nodes per block
    joined = gpd.sjoin(nodes, df, how="inner", predicate="intersects")

    # Group by the census block index and sum the demands
    demand_sum = joined.groupby(joined.index_right)['count'].sum()

    # Merge the aggregated demand back into the census block GeoDataFrame
    df['count'] = df.index.map(demand_sum).fillna(0)
    df['demand_per_node'] = df['peak_LPS']/df['count']

    # Match up nodes to those in df
    new_nodes = gpd.sjoin(nodes, df, how="inner", predicate="intersects")

    # Remove duplicaztes
    combined_nodes = new_nodes[~new_nodes.index.duplicated(keep='first')].copy()

    combined_nodes['demand'] = combined_nodes['demand_per_node'] # Set up demand charecteristic


    # Remove unnesecraty data
    clean_nodes = combined_nodes.drop(columns=['outdoor.demand',"wastewater.supply","indoor_outdoor",
                            "cii.demand","day","count_left","indoor.demand",
                            "year_month","days","LPS","peak_LPS","count_right",
                            "gray_produced","demand_per_node","year","month","id","NAME"])

    BASIC_DEMAND = 0.2 # Fill na values with a chosen demand
    clean_nodes['demand'] = clean_nodes['demand'].fillna(BASIC_DEMAND)

    clean_nodes = clean_nodes.set_index("name") # Set index to non-changing value

    for node in G.nodes: # Copy over data to G
        key = G.nodes[node]['name']
        if key in clean_nodes.index:
            G.nodes[node]['demand'] = clean_nodes['demand'][key]
        else:
            G.nodes[node]['demand'] = BASIC_DEMAND
    
    return G

def attach_IUWM_demands(G,filename,api_key,**kwargs):
    """Runs all functions needed to manipulate a IUWM demand input file and place demands onto a osmnx network

    Args:
        G (osmnx graph): The osmnx graph matching the IUWM data
        filename (str): The location of the IUWM csv
        api_key (str): A Census api key 

    Returns:
        _type_: _description_
    """
    demand_df = read_IUWM(filename)
    # print("\n\n1\n",demand_df.head())
    demand_df = fix_IUWM(demand_df)
    # print("\n\n2\n",demand_df.head())


    demand_df = find_max_month(demand_df)
    # print("\n\n3\n",demand_df.head())

    demand_df = find_demand(demand_df)
    # print("\n\n4\n",demand_df.head())

    census_df = get_census_info(api_key,**kwargs)
    # print("\n\n5\n",census_df.head())


    data_df = match_geoids(demand_df,census_df)
    # print("\n\n6\n",data_df.head())

    data_df = add_geometry_to_census(data_df,**kwargs)
    # print("\n\n7\n",data_df.head())


    G = match_node_to_census(G,data_df)
    return G



