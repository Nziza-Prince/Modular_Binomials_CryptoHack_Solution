from math import gcd

def read_values_from_file(filename):
    """Read the values from a file and return them as variables."""
    values = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                try:
                    # Convert to integer
                    values[key] = int(value)
                except ValueError:
                    # If there's an error, keep as string
                    values[key] = value
    
    return values.get('N'), values.get('e1'), values.get('e2'), values.get('c1'), values.get('c2')

def extended_gcd(a, b):
    """Extended Euclidean Algorithm to find gcd and Bézout coefficients."""
    if a == 0:
        return b, 0, 1
    else:
        gcd, x, y = extended_gcd(b % a, a)
        return gcd, y - (b // a) * x, x

def mod_inverse(a, m):
    """Find the modular multiplicative inverse of a under modulo m."""
    gcd, x, y = extended_gcd(a, m)
    if gcd != 1:
        raise Exception('Modular inverse does not exist')
    else:
        return x % m

def solve_for_primes(N, e1, e2, c1, c2):
    """
    Solve for p and q using the modular binomial method as described.
    
    Given:
    - N = p * q
    - c1 = (2*p + 3*q)^e1 mod N
    - c2 = (5*p + 7*q)^e2 mod N
    
    Where:
    - a1 = 2, b1 = 3
    - a2 = 5, b2 = 7
    """
    print("Solving for p and q using the modular binomial method...")
    
    # For our specific case:
    a1, b1 = 2, 3
    a2, b2 = 5, 7
    
    # Step 1: Calculate modular inverses of a1 and a2
    try:
        a1_inv = mod_inverse(a1, N)
        a2_inv = mod_inverse(a2, N)
    except Exception as e:
        print(f"Error calculating modular inverse: {e}")
        return None, None
    
    # Step 2: Calculate the terms from the solution:
    # q = gcd(pow(a2,(-e2 * e1),N) * pow(c2, e1, N) - pow(a1, (-e1 * e2), N) * pow(c1, e2, N), N)
    
    # Calculate a2^(-e2*e1) mod N
    a2_neg_power = pow(a2_inv, (e2 * e1) % (N-1), N)
    
    # Calculate a1^(-e1*e2) mod N
    a1_neg_power = pow(a1_inv, (e1 * e2) % (N-1), N)
    
    # Calculate c2^e1 mod N
    c2_power = pow(c2, e1, N)
    
    # Calculate c1^e2 mod N
    c1_power = pow(c1, e2, N)
    
    # Calculate the first term: (a2^(-e2*e1) * c2^e1) mod N
    term1 = (a2_neg_power * c2_power) % N
    
    # Calculate the second term: (a1^(-e1*e2) * c1^e2) mod N
    term2 = (a1_neg_power * c1_power) % N
    
    # Calculate the difference and take GCD with N
    q_candidate = gcd((term1 - term2) % N, N)
    
    # If we found a factor, calculate p
    if q_candidate > 1 and N % q_candidate == 0:
        p_candidate = N // q_candidate
        
        # Verify our results
        if p_candidate * q_candidate == N:
            print(f"Found factors:")
            print(f"p = {p_candidate}")
            print(f"q = {q_candidate}")
            
            # Double-check with the original equations
            test_c1 = pow((a1 * p_candidate + b1 * q_candidate), e1, N)
            test_c2 = pow((a2 * p_candidate + b2 * q_candidate), e2, N)
            
            print("\nVerification:")
            print(f"Original c1: {c1}")
            print(f"Computed c1: {test_c1}")
            print(f"Match: {c1 == test_c1}")
            
            print(f"Original c2: {c2}")
            print(f"Computed c2: {test_c2}")
            print(f"Match: {c2 == test_c2}")
            
            return p_candidate, q_candidate
    
    # If the above approach didn't work, try swapping p and q
    p_candidate = q_candidate
    q_candidate = N // p_candidate
    
    if p_candidate * q_candidate == N:
        print(f"Found factors (after swap):")
        print(f"p = {p_candidate}")
        print(f"q = {q_candidate}")
        
        return p_candidate, q_candidate
    
    # Try the same approach but to find p directly
    # p = gcd(pow(b2,(-e2 * e1),N) * pow(c2, e1, N) - pow(b1, (-e1 * e2), N) * pow(c1, e2, N), N)
    
    try:
        b1_inv = mod_inverse(b1, N)
        b2_inv = mod_inverse(b2, N)
        
        b2_neg_power = pow(b2_inv, (e2 * e1) % (N-1), N)
        b1_neg_power = pow(b1_inv, (e1 * e2) % (N-1), N)
        
        term1 = (b2_neg_power * c2_power) % N
        term2 = (b1_neg_power * c1_power) % N
        
        p_candidate = gcd((term1 - term2) % N, N)
        
        if p_candidate > 1 and N % p_candidate == 0:
            q_candidate = N // p_candidate
            
            if p_candidate * q_candidate == N:
                print(f"Found factors (using b coefficients):")
                print(f"p = {p_candidate}")
                print(f"q = {q_candidate}")
                
                return p_candidate, q_candidate
    except Exception as e:
        print(f"Error in alternative approach: {e}")
    
    print("Could not find the factors with the given approach.")
    return None, None

# Main execution
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "data.txt"  # Default filename
        
    print(f"Reading values from file: {filename}")
    
    try:
        N, e1, e2, c1, c2 = read_values_from_file(filename)
        
        if None in (N, e1, e2, c1, c2):
            print("Error: Missing values in the file. Make sure the file contains N, e1, e2, c1, and c2.")
            sys.exit(1)
            
        # Solve for p and q using the modular binomial method
        p, q = solve_for_primes(N, e1, e2, c1, c2)
        
        if p and q:
            # Format the answer as required (this might be for CTF challenges)
            print(f"\nAnswer: p = {p}, q = {q}")
        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)