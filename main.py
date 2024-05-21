import json
import sys
import time

import paramiko

from ProjectionToMap import ProjectionToMap
from projections.BTEDymaxionProjection import BTEDymaxionProjection


def ftp(host: str, port: int, username: str, password: str, world_name: str, output_filename: str):
    transport = paramiko.Transport((host, port))
    transport.connect(username=username, password=password)
    sftp = paramiko.SFTPClient.from_transport(transport)

    with open("data/" + output_filename, "w") as outfile:
        for file in sftp.listdir_iter('./' + world_name + '/region2d'):
            nums = file.filename.split('.')
            x = nums[0]
            y = nums[1]
            if nums[-1] == "2dr":
                outfile.write(f'{x} {y}\n')

    sftp.close()
    transport.close()


def update_chunk_list(regions, host, port, password):
    for region in regions:
        ftp(host, port, region["username"], password, region["world_name"],
            region["region_name"] + "_chunks")


def main(should_update_chunk_list=True, should_prune_chunks=True, should_show_map=False):
    config = None

    print("Loading configuration")

    try:
        with open("config.json") as json_file:
            config = json.load(json_file)
    except IOError:
        print("No configuration file found.")
        sys.exit(1)

    print("Config loaded.")

    if should_update_chunk_list:
        print("Updating chunk list.")
        start_time = time.time()
        update_chunk_list(config["regions"], config["host"], config["port"], config["password"])
        end_time = time.time()
        execution_time = end_time - start_time
        print("Chunk list updated. Execution time:", str(round(execution_time, 2)), "seconds.")
    else:
        print("Skipping chunk list update.")

    print("Loading projection and chunk data.")
    projection = BTEDymaxionProjection()
    print("Projection loaded.")

    map = ProjectionToMap(config["regions"], projection, config["geo"])
    num_chunks = sum([len(region["points"]) for region in map.points_list])
    print(f"Chunk data loaded. Loaded {num_chunks}")

    if should_prune_chunks:
        print("Pruning chunks to region and state boundaries")
        map.points_list = map.get_points_inside()
        num_pruned_chunks = sum([len(region["points"]) for region in map.points_list])
        print(f"Pruned {num_chunks - num_pruned_chunks} chunks.")
    else:
        print("Skipping chunk pruning.")

    if should_show_map:
        print("Opening map")
        map.print()
    else:
        print("Skipping map")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main(False, True, False)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
