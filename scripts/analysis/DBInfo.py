
from pymongo import MongoClient

class DBInfo(dict):
    # Store db connection settings with a bunch of presets. I'm storing the connection settings
    # and not the mongo Collections because you can't pickle Collections
    # Using dot notation access returns mongo Collection objects
    _DB_INFO = {'compil':
                   {'protDB': ("wl-cmadmin.scripps.edu",27018,"ProtDB_072114","ProtDB_072114"),
                    'seqDB': ("wl-cmadmin.scripps.edu",27018,"SeqDB_072114","SeqDB_072114"),
                    'massDB': ("wl-cmadmin.scripps.edu",27018,"MassDB_072114","MassDB_072114"),
                    'taxDB': ("wl-cmadmin.scripps.edu",27017,"TaxDB_072114","TaxDB_072114"),
                    'hashDB': ("wl-cmadmin.scripps.edu",27017,"HashDB_072114","HashDB_072114"),
                    'domainDB': ("wl-cmadmin.scripps.edu",27018,"DomainDB_072114","DomainDB_072114"),
                    'taxDB': ("wl-cmadmin.scripps.edu", 27017,"TaxDB_072114","TaxDB_072114"),
                    'clusterDB': ("wl-cmadmin.scripps.edu", 27017,"071414_ComPIL_forwardonly_0_7","071414_ComPIL_forwardonly_0_7"),
                    'ncbi_taxonomy': ('wl-cmadmin.scripps.edu', 27017, "taxonomy","taxonomy")
                    },
               'compil_mgm':
                   {'protDB': ("wl-cmadmin.scripps.edu",27018,"ProtDB_20151009_compil_mgm","ProtDB_20151009_compil_mgm"),
                    'seqDB': ("wl-cmadmin.scripps.edu",27018,"SeqDB_20151009_compil_mgm","SeqDB_20151009_compil_mgm"),
                    'massDB': ("wl-cmadmin.scripps.edu",27018,"MassDB_20151009_compil_mgm","MassDB_20151009_compil_mgm"),
                    'taxDB': ("wl-cmadmin.scripps.edu",27017,"TaxDB_20151009","TaxDB_20151009"),
                    'hashDB': ("wl-cmadmin.scripps.edu",27017,"HashDB_20151009","HashDB_20151009"),
                    'domainDB': ("wl-cmadmin.scripps.edu",27018,"DomainDB_072114","DomainDB_072114"),
                    'taxDB': ("wl-cmadmin.scripps.edu", 27017,"TaxDB_20151009","TaxDB_20151009"),
                    'clusterDB': ("wl-cmadmin.scripps.edu", 27017,"20151009_compil_mgm_cdhit07","20151009_compil_mgm_cdhit07"),
                    'ncbi_taxonomy': ('wl-cmadmin.scripps.edu', 27017, "taxonomy","taxonomy")
                    },
                # ssh -fNL 27028:localhost:27018 gbwl
                # ssh -fNL 27027:localhost:27017 gbwl
                'compil_mgm_localhost':
                    {'protDB': ("localhost",27028,"ProtDB_20151009_compil_mgm","ProtDB_20151009_compil_mgm"),
                    'seqDB': ("localhost",27028,"SeqDB_20151009_compil_mgm","SeqDB_20151009_compil_mgm"),
                    'massDB': ("localhost",27028,"MassDB_20151009_compil_mgm","MassDB_20151009_compil_mgm"),
                    'taxDB': ("localhost",27027,"TaxDB_20151009","TaxDB_20151009"),
                    'hashDB': ("localhost",27027,"HashDB_20151009","HashDB_20151009"),
                    'domainDB': ("localhost",27028,"DomainDB_072114","DomainDB_072114"),
                    'taxDB': ("localhost", 27027,"TaxDB_20151009","TaxDB_20151009"),
                    'clusterDB': ("localhost", 27027,"20151009_compil_mgm_cdhit07","20151009_compil_mgm_cdhit07"),
                    'ncbi_taxonomy': ('localhost', 27027, "taxonomy","taxonomy")
                    },
                'compil_mgm_food':
                    None
                }

    def __init__(self, preset=None, **kwargs):
        if preset and preset not in DBInfo._DB_INFO:
            raise ValueError("preset must be one of: {}".format(DBInfo._DB_INFO.keys()))
        if preset:
            self.update(DBInfo._DB_INFO[preset])
        # Can give both preset and new args which will override whats in the preset
        self.update(dict(**kwargs))
    
    @staticmethod
    def mongoclient_builder(db_info_tuple):
        return MongoClient(host=db_info_tuple[0], port=db_info_tuple[1])[db_info_tuple[2]][db_info_tuple[3]]

    def to_URI(self, key):
        db_info_tuple = self.__getitem__(key)
        return "mongodb://" + db_info_tuple[0] + ":" + str(db_info_tuple[1])
        
    # overload dot-access to class w/o interfering with __attributes
    def __getattr__(self, key):
        if key.startswith("__"):
            return dict.__getattribute__(self, key)
        db_info_tuple = self[key]
        return MongoClient(host=db_info_tuple[0], port=db_info_tuple[1])[db_info_tuple[2]][db_info_tuple[3]]
        #return self.__class__.mongoclient_builder(self[key])

        
if __name__ == "__main__":
    # from metaproteomics.analysis.DBInfo import DBInfo
    import dill
    db_info = DBInfo('compil', protDB = ("wl-cmadmin.scripps.edu",27018,"ProtDB_new","ProtDB_new"))
    print(db_info)
    print(db_info.protDB)
    print(db_info.protDB.database.name)
    with open('test123','wb') as f:
        dill.dump(db_info, f)

""" # w/o overloading __getattr__ ?
    @staticmethod
    def mongoclient_builder(db_info_tuple):
        return MongoClient(host=db_info_tuple[0], port=db_info_tuple[1])[db_info_tuple[2]][db_info_tuple[3]]
    
    @staticmethod
    def mongoURI_builder(db_info_tuple):
        return "mongodb://" + db_info_tuple[0] + ":" + str(db_info_tuple[1])

    def to_URI(self, key):
        return DBInfo.mongoURI_builder(self.__getitem__(key))
    
    def to_mongo(self, key):
        return DBInfo.mongoclient_builder(self.__getitem__(key))
"""