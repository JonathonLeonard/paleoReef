import pandas as pd
import pygplates
import sys

class PaleoRotator:
    def __init__(self, rotation_file, static_polygon_file, japan_plate_ids=None, honshu_plate_id=629):
        """
        Initializes the PaleoRotator with the path to the rotation file.

        :param rotation_file: Path to the GPlates rotation file.
        """
        self.rotation_file = rotation_file
        self.static_polygon_file = static_polygon_file
        self.rotation_model = pygplates.RotationModel(self.rotation_file)
        self.japan_plate_ids = japan_plate_ids or list(range(619, 632))
        self.honshu_plate_id = honshu_plate_id

    def moveJapan2Honshu(self, gpml_feature_collection):
        for feature in gpml_feature_collection:
            if feature.get_reconstruction_plate_id() in self.japan_plate_ids:
                feature.set_reconstruction_plate_id(self.honshu_plate_id)
        return gpml_feature_collection
    
    def points2gpml(self, lons, lats, maxAges=None, minAges=None, names=None):
        input_points = [pygplates.PointOnSphere(lat,lon) for lon, lat in zip(lons, lats)]
        point_features = []
        for i, point in enumerate(input_points):
            point_feature = pygplates.Feature(pygplates.FeatureType.gpml_unclassified_feature)
            point_feature.set_geometry(point)
            if maxAges is not None and minAges is not None:
                print(f"Setting valid time for point {i}: maxAge={maxAges[i]}, minAge={minAges[i]}")
                point_feature.set_valid_time(maxAges[i], minAges[i])
            if names is not None:
                point_feature.set_name(str(names[i]))
            point_features.append(point_feature)

        properties2copy = [pygplates.PartitionProperty.reconstruction_plate_id]
        # if maxAges is not None and minAges is not None:
        #     properties2copy.append(pygplates.PartitionProperty.valid_time_period)
        
        assigned_point_features = pygplates.partition_into_plates(
            self.static_polygon_file,
            self.rotation_model,
            point_features,
            properties_to_copy = properties2copy,
            reconstruction_time = 0
        )
        assigned_point_feature_collection = pygplates.FeatureCollection(assigned_point_features)
        assigned_point_feature_collection = self.moveJapan2Honshu(assigned_point_feature_collection)
        return assigned_point_feature_collection
    
    def rotatePoints(self, gpml_feature_collection, rotation_times, anchor_plate_id=0):
        reconstructed_lats = []
        reconstructed_lons = []
        
        for i, pt in enumerate(gpml_feature_collection):
            age = rotation_times[i]
            reconstructed_pts = []
            print(self.rotation_model)
            print(f'Reconstructing point {i} at age {age} with anchor plate ID {anchor_plate_id}')
            # Reconstruct the point at the specified age
            pygplates.reconstruct(pt, self.rotation_model, reconstructed_pts, reconstruction_time=age, anchor_plate_id=anchor_plate_id)
            try:
                latlons = reconstructed_pts[0].get_reconstructed_geometry().to_lat_lon_list()
                # Print detailed information about the point and reconstruction results for debugging
                print(f"Plate ID: {pt.get_reconstruction_plate_id()}")
                print(f"Valid time: {pt.get_valid_time()}")
                print(f"Original geometry (lat, lon): {pt.get_geometry().to_lat_lon()}")
                print(f"Reconstruction time: {age}")
                print(f"Reconstructed geometry (lat, lon): {latlons[0]}")
                print(f"Reconstructed points: {reconstructed_pts}")
                reconstructed_lats.append(round(latlons[0][0],2))
                reconstructed_lons.append(round(latlons[0][1],2))
            except Exception as e:
                print(e)
                print(f'Error reconstructing point {i} at age {age}. Skipping.')
                # Print the reconstructed points for debugging
                # Print detailed information about the point and reconstruction results for debugging
                print(f"Plate ID: {pt.get_reconstruction_plate_id()}")
                print(f"Valid time: {pt.get_valid_time()}")
                print(f"Original geometry (lat, lon): {pt.get_geometry().to_lat_lon()}")
                print(f"Reconstructed points: {reconstructed_pts}")
                sys.exit()
                reconstructed_lats.append(None)
                reconstructed_lons.append(None)
            print('\n')
        return reconstructed_lats, reconstructed_lons

