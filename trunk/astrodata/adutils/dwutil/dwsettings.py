from astrodata import Lookups

package_classes = Lookups.compose_multi_table(
                            "*/warehouse_settings", "warehouse_package")
package_dict = {}
if "warehouse_package" in package_classes:
    warehouse_packages = package_classes["warehouse_package"]
    for package in warehouse_packages:
        for key in package:
            if key not in package_dict:
                package_dict[key] = package[key]
else:
    warehouse_packages = []

dataset_extensions_dict = Lookups.compose_multi_table(
                        "*/filetypes", "data_object_precedence")

if "data_object_precedence" in dataset_extensions_dict:
    dataset_extensions = dataset_extensions_dict["data_object_precedence"]
else:
    dataset_extensions = None

ingest_sources_dict = Lookups.compose_multi_table(
                            "*/warehouse_daemon.py",
                            "ingest_sources")
if "ingest_sources" in ingest_sources_dict:
    ingest_sources = ingest_sources_dict["ingest_sources"]
else:
    ingest_sources = None
    