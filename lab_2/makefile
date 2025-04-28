# Определяем компилятор и флаги 
CXX = g++ 
CXXFLAGS = -Wall -Wextra -std=c++11 

# Определяем целевую программу 
TARGET = main 

# Определяем исходные файлы и объектные файлы 
SRCS = 2.cpp
OBJS = $(SRCS:.cpp=.o) 

# Правило по умолчанию 
all: $(TARGET) 

# Правило для создания исполняемого файла 
$(TARGET): $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $^ 

# Правило для компиляции объектных файлов 
%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

# Правило для очистки 
clean: 
	del /Q $(OBJS) $(TARGET)

.PHONY: all clean
