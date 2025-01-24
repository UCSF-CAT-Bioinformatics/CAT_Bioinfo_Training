# Course Title: Bash Scripting in a Flash (with AI Power!)

A 1-hour introductory Bash scripting course, culminating in using AI assistance:

**Course Objectives:**

*   Understand what Bash scripting is and why it's useful.
*   Learn fundamental Bash commands and syntax.
*   Write simple Bash scripts.
*   Discover how AI tools like Gemini can aid in script development.

## Introduction & What is Bash? 

A command-line interpreter and scripting language for Unix-like systems (Linux, macOS).

Highlighting its uses: automating tasks, system administration, file manipulation, etc.

Knowing how to use Bash and script in bash provides give you power over your system and continues to be relevant even in today's GUI-driven world.

## Recalling Basic Bash Commands & Syntax

*   **Navigation:**
    `pwd` (print working directory), `cd` (change directory), `ls` (list files).
*   **File Manipulation:**
    `touch` (create file), `mkdir` (make directory), `rm` (remove), `cp` (copy), `mv` (move).
*   **Viewing Files:**
    `cat` (concatenate and display), `less` (page through), `head` (first few lines), `tail` (last few lines).
*   **Input/Output Redirection:**
    `>` (output to file, overwriting), `>>` (append to file), `<` (input from file), `|` (pipe output to another command).

**Example:**
*   Create a directory, create a file inside it, write some text to the file (echo, nano), and display its contents.

Within a bash script you can include anything that can be used on the command line ('bash')

## Writing a Simple Bash Script

*   **Shebang:**
    The `#!/bin/bash` line at the beginning of a script (tells the system which interpreter to use).
*   **Variables:**
    Declaring and using variables (e.g., `name="John"`, `echo "Hello, $name"`).
*   **Command lines parameters:**
    *   $@: this contains all the input arguments
    *   $#: the number of arguments passed to the script
    *   $0: the name of the script itself
    *   $1: the first parameter
    *   $2: the second parameter, etc.
*   **Echo:**
    The `echo` command for displaying output.
*   **Comments:**
    Comments are added with `#`.
*   **Permissions:**
    Explain `chmod +x script.sh` to make the script executable.

### **Simple Script Example:**

```bash
#!/bin/bash
name="User"
echo "Hello, $name!"
echo "Today is $(date)"
```

*   **Hands-on exercise:** Create a simple script that prints your name and the current date, BUT your name is provide on the command line.

## Using AI (Gemini) for Bash Scripting, AI Assistance in Coding

So we've learned the basics of Bash scripting. Now, let's talk about a powerful tool that can make our coding lives much easier: Artificial Intelligence, specifically models like Gemini, or ChatGPT.

### What is Gemini?

Gemini is a large language model, which means it's been trained on a massive amount of text and code. This allows it to understand human language and, importantly, programming languages like Bash.

### How Gemini Helps with Coding - Code Generation

One of the most useful things Gemini can do is generate code for you. Let's say you need a script to perform a specific task, but you're not quite sure how to write it. You can simply describe what you want in plain English, and Gemini can often generate the Bash code for you.

**(Demonstration 1: Simple Task)**

Let's try a simple example. I'll ask Gemini: "Write a Bash script to list all files in the current directory that are larger than 1MB."

_Gemini generates the code. Example output might be similar to this, but will vary:_

```bash
#!/bin/bash

find . -type f -size +1M -print
```

Okay, it gave us this code. Let's break it down. `find` is a command for searching files. `-type f` means we're looking for files. `-size +1M` means files larger than 1 megabyte. And `-print` tells it to display the results. Pretty neat, right?

### How Gemini Helps with Coding - Code Explanation

Gemini isn't just about generating code. It can also explain existing code. If you come across a Bash script and you're not sure what it does, you can give it to Gemini, and it can provide an explanation.

**Demonstration 2: Code Explanation**

Let's take a slightly more complex example. Let's say we have this script:

_Lets take the following script:_

```bash
#!/bin/bash

for file in *; do
  if [[ -f "$file" ]]; then
    size=$(wc -c < "$file")
    echo "$file: $size bytes"
  fi
done
```

Now lets ask Gemini: "Explain this Bash script:" (paste the script).

**Gemini provides an explanation, something like:**

"This script iterates through all files in the current directory. For each file, it checks if it's a regular file (not a directory). If it is, it calculates the file size in bytes using `wc -c` and then prints the filename and its size."

## How Gemini Helps with Coding - Code Improvement/Refactoring

Another powerful feature is code improvement or refactoring. You can provide a script and ask Gemini to make it more efficient, readable, or to fix potential bugs.

**Demonstration 3: Code Improvement - Example**

Let's say we have the previous script, and we want to make it handle spaces in filenames correctly. We could ask Gemini: "Improve this Bash script to handle filenames with spaces:" (paste the script). 

Gemini might suggest quoting the `$file` variable in the `wc -c` command.

## Prompt Engineering

The key to getting good results from Gemini is *prompt engineering*. This means crafting your prompts effectively.

*   **Be specific:** The more detail you give, the better the results.
*   **Provide context:** If the task is related to a specific project, provide some background.
*   **Iterate:** If the first result isn't perfect, rephrase your prompt and try again.

## Important Considerations

It's essential to understand that Gemini isn't a replacement for learning to code. It's a tool to assist you.

*   **Review and Understand:** Always review the code Gemini generates. Make sure you understand how it works.
*   **Testing:** Always test the generated code thoroughly.
*   **Not Always Perfect:** Gemini may sometimes produce incorrect or suboptimal code.

## Real Example within the CAT

Lets create a run characteristic report!

**Materials:**

*   Text editor (e.g., Nano, Vim, VS Code) for demonstrations.
*   Access to a terminal (Linux/macOS or Windows Subsystem for Linux (WSL)).
*   Internet access for demonstrating Gemini.
