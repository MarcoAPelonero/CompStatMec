import sys

def main():
    print("Python script is ready to receive input.")
    for line in sys.stdin:
        print(f"Received: {line.strip()}")
        print("Processed: " + line.strip()[::-1])  # Example processing: reversing the input

if __name__ == "__main__":
    main()
