def append_stl_file(source_file, destination_file):
    try:
        with open(source_file, 'rb') as source:
            source_data = source.read()
        
        with open(destination_file, 'ab') as destination:
            destination.write(source_data)
        
        print(f"Appended {source_file} to {destination_file} successfully.")
    
    except FileNotFoundError:
        print("One or both of the files does not exist.")
    
    except Exception as e:
        print(f"An error occurred: {e}")

append_stl_file('flat_plate_2og.stl', 'flat_plate_1.stl')